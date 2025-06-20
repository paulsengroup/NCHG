// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: GPL-3.0
//
// This library is free software: you can redistribute it and/or
// modify it under the terms of the GNU Public License as published
// by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Public License along
// with this library.  If not, see
// <https://www.gnu.org/licenses/>.

#include "nchg/parquet_stats_file_merger.hpp"

#include <arrow/ipc/writer.h>
#include <arrow/memory_pool.h>
#include <arrow/type.h>
#include <arrow/util/base64.h>
#include <duckdb.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/std.h>
#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <hictk/chromosome.hpp>
#include <memory>
#include <numeric>
#include <random>
#include <ranges>
#include <vector>

#include "nchg/metadata.hpp"
#include "nchg/parquet_stats_file_reader.hpp"
#include "nchg/parquet_stats_file_writer.hpp"
#include "nchg/version.hpp"

// NOLINTBEGIN(*-dcl58-cpp, *-owning-memory)
template <>
struct std::default_delete<duckdb_config> {
  constexpr void operator()(duckdb_config* cfg) const noexcept {
    duckdb_destroy_config(cfg);
    delete cfg;
  }
};

template <>
struct std::default_delete<duckdb_database> {
  constexpr void operator()(duckdb_database* db) const noexcept {
    duckdb_close(db);
    delete db;
  }
};

template <>
struct std::default_delete<duckdb_connection> {
  constexpr void operator()(duckdb_connection* con) const noexcept {
    duckdb_disconnect(con);
    delete con;
  }
};

template <>
struct std::default_delete<duckdb_result> {
  constexpr void operator()(duckdb_result* res) const noexcept {
    duckdb_destroy_result(res);
    delete res;
  }
};

template <>
struct std::default_delete<duckdb_prepared_statement> {
  constexpr void operator()(duckdb_prepared_statement* stmt) const noexcept {
    duckdb_destroy_prepare(stmt);
    delete stmt;
  }
};
// NOLINTEND(*-dcl58-cpp, *-owning-memory)

namespace nchg {

[[nodiscard]] static std::vector<std::vector<std::filesystem::path>> preprocess_file_list(
    const phmap::btree_map<hictk::Chromosome, std::vector<std::filesystem::path>>& files) {
  using Buffer = std::vector<std::vector<std::filesystem::path>>;
#ifdef __cpp_lib_ranges_to_container
  return files | std::ranges::views::values | std::ranges::to<Buffer>();
#else
  const auto values = files | std::ranges::views::values;
  return Buffer{values.begin(), values.end()};
#endif
}

class DuckDBHandle {
  std::unique_ptr<duckdb_database> _db;
  std::unique_ptr<duckdb_connection> _con;

 public:
  DuckDBHandle(const std::filesystem::path& tmpdir, std::uint64_t threads,
               std::uint64_t memory_limit_mb)
      : _db(open_database(make_config(tmpdir, memory_limit_mb, threads))),
        _con(open_connection(*_db)) {}
  DuckDBHandle(const DuckDBHandle&) = delete;
  DuckDBHandle(DuckDBHandle&&) noexcept = delete;

  ~DuckDBHandle() noexcept = default;

  DuckDBHandle& operator=(const DuckDBHandle&) = delete;
  DuckDBHandle& operator=(DuckDBHandle&&) noexcept = delete;

  [[nodiscard]] const duckdb_database& db() const noexcept {
    assert(_db);
    return *_db;
  }
  [[nodiscard]] duckdb_database& db() noexcept {
    assert(_db);
    return *_db;
  }

  [[nodiscard]] const duckdb_connection& connection() const noexcept {
    assert(_con);
    return *_con;
  }
  [[nodiscard]] duckdb_connection& connection() noexcept {
    assert(_con);
    return *_con;
  }

  [[nodiscard]] std::unique_ptr<duckdb_prepared_statement> prepare_statement(
      const std::string& str) {
    SPDLOG_DEBUG("DuckDB: about to prepare the following statement:\n{}", str);
    auto stmt = std::make_unique_for_overwrite<duckdb_prepared_statement>();
    if (duckdb_prepare(connection(), str.c_str(), stmt.get()) == DuckDBError) {
      throw std::runtime_error(
          fmt::format("failed to prepare statement: {}", duckdb_prepare_error(*stmt)));
    }
    return stmt;
  }

  void query(duckdb_prepared_statement& stmt) {
    auto res = std::make_unique_for_overwrite<duckdb_result>();
    if (duckdb_execute_prepared(stmt, res.get()) != DuckDBSuccess) {
      // calling duckdb_result_error in this way is ok because the message is stored inside the res
      // object
      throw std::runtime_error(
          fmt::format("failed to execute SQL query: {}", duckdb_result_error(res.get())));
    }
  }

  void query(const std::string& stmt) {
    auto res = std::make_unique_for_overwrite<duckdb_result>();
    if (duckdb_query(*_con, stmt.c_str(), res.get()) != DuckDBSuccess) {
      // calling duckdb_result_error in this way is ok because the message is stored inside the res
      // object
      throw std::runtime_error(
          fmt::format("failed to execute SQL query: {}", duckdb_result_error(res.get())));
    }
  }

 private:
  [[nodiscard]] static std::unique_ptr<duckdb_config> make_config(
      const std::filesystem::path& tmpdir, std::uint64_t memory_limit_mb,
      std::uint64_t num_threads) {
    assert(num_threads > 0);
    try {
      auto cfg = std::make_unique_for_overwrite<duckdb_config>();
      if (duckdb_create_config(cfg.get()) != DuckDBSuccess) {
        throw std::runtime_error("failed to allocate memory for the config object");
      }
      assert(std::filesystem::is_directory(tmpdir));

      auto set_or_throw = [&](const auto& key, const auto& value) {
        const auto str_value = fmt::to_string(value);
        if (duckdb_set_config(*cfg, key, str_value.c_str()) != DuckDBSuccess) {
          throw std::runtime_error(fmt::format(R"(failed to set "{}" to "{}")", key, str_value));
        }
      };

      set_or_throw("temp_directory", tmpdir.string());
      set_or_throw("memory_limit", fmt::format("{}MB", memory_limit_mb));
      set_or_throw("preserve_insertion_order", true);
      set_or_throw("allow_community_extensions", false);
      set_or_throw("threads", num_threads);
      set_or_throw("allocator_background_threads", true);

      return cfg;
    } catch (const std::exception& e) {
      throw std::runtime_error(fmt::format("failed to initialize DucDB's config: {}", e.what()));
    }
  }

  [[nodiscard]] static std::unique_ptr<duckdb_database> open_database(
      const std::unique_ptr<duckdb_config>& cfg) {
    try {
      char* error_msg{};
      auto db = std::make_unique_for_overwrite<duckdb_database>();
      if (duckdb_open_ext(nullptr, db.get(), *cfg, &error_msg) != DuckDBSuccess) {
        if (error_msg) {
          std::runtime_error e{error_msg};
          duckdb_free(error_msg);
          throw e;
        }
        throw std::runtime_error("unknown error");
      }

      return db;
    } catch (const std::exception& e) {
      throw std::runtime_error(fmt::format("failed to initialize DuckDB's database: {}", e.what()));
    }
  }

  [[nodiscard]] static std::unique_ptr<duckdb_connection> open_connection(duckdb_database& db) {
    auto con = std::make_unique_for_overwrite<duckdb_connection>();
    if (duckdb_connect(db, con.get()) != DuckDBSuccess) {
      throw std::runtime_error("failed to connect to DuckDB database");
    }
    return con;
  }
};

[[nodiscard]] static std::string try_read_arrow_schema(
    const std::vector<std::vector<std::filesystem::path>>& files) {
  std::shared_ptr<arrow::Schema> schema;
  for (const auto& file : files | std::ranges::views::join) {
    if (!schema) {
      schema = ParquetStatsFileReader::get_file_schema(file);
      continue;
    }

    const auto new_schema = ParquetStatsFileReader::get_file_schema(file);
    if (new_schema && *schema != *new_schema) {
      throw std::runtime_error("caught an attempt to merge files with different schemas");
    }
  }

  if (!schema) {
    return "";
  }

  // Based on
  // https://github.com/apache/arrow/blob/79281cc7e815dc10c1cfcebfd900c12924ddadf2/cpp/src/parquet/arrow/writer.cc#L526-L551
  const auto res = arrow::ipc::SerializeSchema(*schema, arrow::default_memory_pool());
  if (!res.ok()) {
    return "";
  }

  const auto schema_as_string = res.ValueUnsafe()->ToString();
  return arrow::util::base64_encode(schema_as_string);
}

[[nodiscard]] static hictk::Reference read_chromosomes_checked(
    const std::vector<std::vector<std::filesystem::path>>& files) {
  hictk::Reference chroms;
  for (const auto& f : files | std::ranges::views::join) {
    if (chroms.empty()) {
      chroms = *ParquetStatsFileReader{f, ParquetStatsFileReader::RecordType::infer}.chromosomes();
      continue;
    }

    if (chroms !=
        *ParquetStatsFileReader{f, ParquetStatsFileReader::RecordType::infer}.chromosomes()) {
      throw std::runtime_error("caught attempt to merge files using different reference genomes");
    }
  }
  return chroms;
}

[[nodiscard]] static std::unique_ptr<duckdb_prepared_statement> declare_chromosome_enum_sql(
    DuckDBHandle& db, const hictk::Reference& chroms) {
  if (chroms.empty()) {
    return {};
  }
  constexpr auto* fmt_string = R"(
CREATE TYPE chromosomes AS ENUM (
  '{}'
);
)";

  return db.prepare_statement(fmt::format(
      fmt_string,
      fmt::join(chroms.remove_ALL() |
                    std::ranges::views::transform([](const auto& chrom) { return chrom.name(); }),
                "',\n  '")));
}

// It's important that this function is pre-declared on top of all the other overloads
static void escape_json_field_for_sql(glz::json_t::val_t& field);

static void escape_json_field_for_sql(glz::json_t& field) { escape_json_field_for_sql(field.data); }

static void escape_json_field_for_sql(glz::json_t::array_t& field) {
  for (auto& f : field) {
    escape_json_field_for_sql(f.data);
  }
}

static void escape_json_field_for_sql(glz::json_t::object_t& field) {
  for (auto& f : field | std::ranges::views::values) {
    escape_json_field_for_sql(f.data);
  }
}

static constexpr void escape_json_field_for_sql([[maybe_unused]] bool field) {}

static constexpr void escape_json_field_for_sql([[maybe_unused]] double field) {}

static void escape_json_field_for_sql(std::string& s) {
  constexpr char q = '\'';
  const auto num_single_quotes = static_cast<std::size_t>(std::ranges::count(s, q));
  std::string buff;
  buff.reserve(s.size() + num_single_quotes);
  if (num_single_quotes == 0) {
    buff.append(s);
  } else {
    for (const auto c : s) {
      if (c == q) {
        buff.append("''");
      } else {
        buff.append(1, c);
      }
    }
  }
  std::swap(s, buff);
}

static void escape_json_field_for_sql(glz::json_t::val_t& field) {
  std::visit(
      [&]<typename T>(T& field_) {
        if constexpr (std::is_same_v<glz::json_t::null_t, std::remove_cvref_t<T>>) {
          field = "NULL";
        } else {
          escape_json_field_for_sql(field_);
        }
      },
      field);
}

[[nodiscard]] static std::string escape_metadata_json_string_sql(std::string_view metadata) {
  if (metadata.empty()) {
    return "{}";
  }

  auto json = parse_json_string(metadata);
  escape_json_field_for_sql(json);
  constexpr glz::opts opts{
      .format = glz::JSON,
      .prettify = true,
      .indentation_width = 6,
      .new_lines_in_arrays = true,
      .raw_string = true,
  };

  std::string buff;
  if (const auto ec = glz::write<opts>(json, buff); ec) {
    throw std::runtime_error(glz::format_error(ec));
  }

  if (buff.empty()) {
    return "{}";
  }

  buff.resize(buff.size() - 1);
  buff.append("    }");

  return buff;
}

[[nodiscard]] static std::unique_ptr<duckdb_prepared_statement> merge_parquets_sql(
    DuckDBHandle& db, const std::vector<std::vector<std::filesystem::path>>& files,
    const std::filesystem::path& dest, std::string_view metadata,
    std::string_view compression_method, std::uint8_t compression_level) {
  // clang-format off
  constexpr auto* statement_template = R"(
COPY (
  SELECT
    * REPLACE (
      -- cast to ENUM to ensure dictionary contains values for all chromosomes
      chrom1::chromosomes AS chrom1,
      chrom2::chromosomes AS chrom2
    )
  FROM
    read_parquet(
      [
        ${}
      ]
    )
)
TO '{}' (
  FORMAT PARQUET,
  COMPRESSION {},
  COMPRESSION_LEVEL {},
  KV_METADATA
    {}
);
)";
  // clang-format on

  const auto num_parquet_files = std::accumulate(
      files.begin(), files.end(), 0UZ, [&](auto acc, const auto& v) { return acc + v.size(); });

  const auto stmt = fmt::format(
      // clang-format off
    statement_template,
    fmt::join(std::ranges::views::iota(1UZ, num_parquet_files + 1), ",\n        $"),
    dest.string(),
    compression_method,
    compression_level,
    metadata
      // clang-format on
  );

  auto pstmt = db.prepare_statement(stmt);

  std::size_t idx = 1;
  for (const auto& file : files | std::ranges::views::join) {
    duckdb_bind_varchar(*pstmt, idx++, file.c_str());
  }

  return pstmt;
}

[[nodiscard]] static std::unique_ptr<duckdb_prepared_statement> sort_and_merge_parquets_sql(
    DuckDBHandle& db, const std::vector<std::vector<std::filesystem::path>>& files,
    const std::filesystem::path& dest, std::string_view metadata,
    std::string_view compression_method, std::uint8_t compression_level) {
  // TODO compression tunable

  std::size_t i = 1;
  std::vector<std::string> partial_queries(files.size());
  // clang-format off
  constexpr auto* statement1_template = R"(
    (
      SELECT
        * REPLACE (
          -- cast to ENUM to ensure dictionary contains values for all chromosomes
          chrom1::chromosomes AS chrom1,
          chrom2::chromosomes AS chrom2
        )
        FROM
          read_parquet(
            [
              ${}
            ]
          )
        -- sort by genomic coordinates
        ORDER BY chrom1, start1, end1, chrom2, start2, end2
    )
)";

  constexpr auto* statement2_template = R"(
COPY (
  SELECT * FROM (
    {}
  )
)
TO '{}' (
  FORMAT PARQUET,
  COMPRESSION {},
  COMPRESSION_LEVEL {},
  KV_METADATA
    {}
);
)";
  // clang-format on

  std::ranges::transform(
      files, partial_queries.begin(), [&i, statement1_template](const auto& files) {
        const auto i0 = i;
        const auto i1 = i0 + files.size();
        i = i1;

        return fmt::format(statement1_template,
                           fmt::join(std::ranges::views::iota(i0, i1), ",\n              $"));
      });

  const auto stmt = fmt::format(statement2_template, fmt::join(partial_queries, "    UNION ALL"),
                                dest.string(), compression_method, compression_level, metadata);

  auto pstmt = db.prepare_statement(stmt);

  i = 1;
  for (const auto& file : files | std::ranges::views::join) {
    duckdb_bind_varchar(*pstmt, i++, file.c_str());
  }

  return pstmt;
}

[[nodiscard]] static std::string generate_metadata_str(
    const hictk::Reference& chroms, const std::vector<glz::json_t>& input_metadata) {
  glz::json_t metadata{
      {"chromosomes", parse_json_string(to_json_string(chroms.remove_ALL()))},
      {"command", "merge"},
      {"date", fmt::format("{:%FT%T}", fmt::gmtime(std::chrono::system_clock::now()))},
      {"input-metadata", parse_json_string(to_json_string(input_metadata))},
      {"version", config::version::str()}};
  return to_json_string(metadata);
}

[[nodiscard]] static std::string generate_metadata_str(
    const std::vector<std::vector<std::filesystem::path>>& files, const hictk::Reference& chroms,
    std::string_view arrow_schema) {
  std::vector<glz::json_t> old_metadata;
  std::ranges::transform(
      files | std::ranges::views::join, std::back_inserter(old_metadata), [&](const auto& path) {
        try {
          const auto metadata = ParquetStatsFileReader::read_metadata(path, {"chromosomes"});
          return parse_json_string(metadata);
        } catch (const std::exception& e) {
          throw std::runtime_error(
              fmt::format("failed to read NCHG metadata from file \"{}\": {}", path, e.what()));
        } catch (...) {
          throw std::runtime_error(
              fmt::format("failed to read NCHG metadata from file \"{}\": unknown error", path));
        }
      });

  // TODO this should be compressed when appropriate
  const auto metadata =
      escape_metadata_json_string_sql(generate_metadata_str(chroms, old_metadata));

  constexpr auto* fmt_string_wo_schema =
      R"({{
'NCHG:format-version': {},
'NCHG:metadata': '{}',
'NCHG:metadata-compression': '{}',
'NCHG:metadata-size': {}
}})";

  constexpr auto* fmt_string_with_schema =
      R"({{
'NCHG:format-version': {},
'NCHG:metadata': '{}',
'NCHG:metadata-compression': '{}',
'NCHG:metadata-size': {},
'ARROW:schema': '{}'
}})";

  if (arrow_schema.empty()) {
    return fmt::format(fmt_string_wo_schema, ParquetStatsFileMerger::format_version, metadata,
                       "None", metadata.size());
  }
  return fmt::format(fmt_string_with_schema, ParquetStatsFileMerger::format_version, metadata,
                     "None", metadata.size(), arrow_schema);
}

ParquetStatsFileMerger::ParquetStatsFileMerger(
    const phmap::btree_map<hictk::Chromosome, std::vector<std::filesystem::path>>& files)
    : _file_groups(preprocess_file_list(files)) {}

void ParquetStatsFileMerger::merge(const std::filesystem::path& dest, bool sort,
                                   const Params& params) {
  if (std::filesystem::exists(dest)) {
    throw std::runtime_error("file already exists");
  }

  DuckDBHandle db{params.tmpdir, params.threads, params.memory_limit_mb};

  // TODO can these two and the metadata scan be consolidated?
  const auto arrow_schema = try_read_arrow_schema(_file_groups);
  const auto chroms = read_chromosomes_checked(_file_groups);

  const std::string metadata =
      arrow_schema.empty() ? "" : generate_metadata_str(_file_groups, chroms, arrow_schema);

  auto create_enum_stmt = declare_chromosome_enum_sql(db, chroms);
  auto merge_parquets_stmt =
      sort ? sort_and_merge_parquets_sql(db, _file_groups, dest, metadata,
                                         params.compression_method, params.compression_level)
           : merge_parquets_sql(db, _file_groups, dest, metadata, params.compression_method,
                                params.compression_level);

  db.query(*create_enum_stmt);
  db.query(*merge_parquets_stmt);
}

}  // namespace nchg
