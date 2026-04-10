#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <variant>

#include "src/ps_full_intersection.cpp"
#include "src/ps_threshold_union.cpp"
#include "src/ps_utils.cpp"

using namespace fulgor;

template <typename FulgorIndex>
struct named_fastq_query_reader {
    named_fastq_query_reader(std::string& query_filename, uint64_t num_threads, FulgorIndex&)
        : rparser({query_filename}, num_threads, num_threads - 1) {
        rparser.start();
    }

    struct query_t {
        std::string name;
        std::string seq;
    };

    struct query_group {
        explicit query_group(named_fastq_query_reader* reader_)
            : reader(reader_), rg(reader->rparser.getReadGroup()) {}

        bool has_next() { return curr_record != rg.end(); }

        void next() { ++curr_record; }

        void value(query_t& query) {
            query.name = curr_record->name;
            query.seq = curr_record->seq;
        }

        bool refill() {
            const bool result = reader->rparser.refill(rg);
            if (result) curr_record = rg.begin();
            return result;
        }

    private:
        named_fastq_query_reader* reader;
        fastx_parser::ReadGroup<klibpp::KSeq> rg;
        std::vector<klibpp::KSeq>::iterator curr_record;
    };

    query_group get_query_group() { return query_group(this); }

    ~named_fastq_query_reader() { rparser.stop(); }

private:
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser;
};

template <typename FulgorIndex>
struct cobs_output_buffer {
    explicit cobs_output_buffer(FulgorIndex const& index_, std::ofstream& file_, std::mutex& mut_)
        : index(index_), file(file_), mut(mut_) {}

    void write(std::string const& query_name, std::vector<pseudoalignment_match> const& matches) {
        buffer << '*' << query_name << '\t' << matches.size() << '\n';
        for (auto const& match : matches) {
            auto stem = std::filesystem::path(std::string(index.filename(match.color))).stem().string();
            buffer << '_' << stem << '\t' << match.score << '\n';
        }
        if (buffer.tellp() > static_cast<std::streamoff>(1 << 14)) flush();
    }

    ~cobs_output_buffer() { flush(); }

private:
    void flush() {
        const std::string outs = buffer.str();
        if (outs.empty()) return;
        std::lock_guard<std::mutex> lock(mut);
        file.write(outs.data(), outs.size());
        buffer.str("");
        buffer.clear();
    }

    FulgorIndex const& index;
    std::ofstream& file;
    std::mutex& mut;
    std::stringstream buffer;
};

template <typename FulgorIndex, typename Formatter, typename QueryReader>
int pseudoalign_worker(FulgorIndex const& index, QueryReader& query_reader,
                       Formatter& formatter, const double threshold, ps_options& options)  //
{
    auto output_buffer = formatter.buffer();
    std::vector<uint32_t> tmp, colors;  // result of pseudoalignment
    std::vector<uint32_t> color_set_ids;
    std::stringstream ss;

    auto qg = query_reader.get_query_group();
    while (qg.refill()) {
        while (qg.has_next()){
            typename QueryReader::query_t query;
            qg.value(query);

            switch (options.algo) {
                case pseudoalignment_algorithm::FULL_INTERSECTION:
                    index.pseudoalign_full_intersection(query.cids, colors, tmp);
                    break;
                case pseudoalignment_algorithm::THRESHOLD_UNION:
                    index.pseudoalign_threshold_union(query.seq, colors, threshold);
                    break;
                default:
                    break;
            }

            if constexpr (std::is_same_v<preprocessed_query_reader, QueryReader>) {
                options.increment_processed_reads(query.ids.size());
                for (auto qid: query.ids) {
                    output_buffer.write(qid, colors);
                    if (!colors.empty()) {
                        options.increment_mapped_reads();
                    }
                }
            } else {
                options.increment_processed_reads();
                output_buffer.write(query.id, colors);
                if (!colors.empty()) {
                    options.increment_mapped_reads();
                }
            }

            colors.clear();
            qg.next();
        }
    }
    return 0;
}

template <typename FulgorIndex, typename Formatter, typename QueryReader>
int pseudoalign_orchestrator(FulgorIndex& index, QueryReader& query_reader,
                Formatter& formatter, const double threshold, ps_options& options) {
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> t;
    t.start();

    uint64_t num_threads = options.num_threads;
    assert(num_threads >= 2);

    if (options.verbose) essentials::logger("*** START: pseudoalignment");
    std::vector<std::thread> workers;
    workers.reserve(num_threads);
    for (uint64_t i = 1; i != num_threads; ++i) {
        workers.push_back(std::thread([&index, &query_reader, &formatter, threshold, &options]() {
            pseudoalign_worker(index, query_reader, formatter, threshold, options);
        }));
    }

    for (auto& w : workers) w.join();

    t.stop();
    if (options.verbose) essentials::logger("*** DONE: pseudoalignment");

    if (options.verbose) {
        std::cout << "processed " << options.num_reads << " reads" << std::endl;
        std::cout << "elapsed = " << t.elapsed() << " millisec / ";
        std::cout << t.elapsed() / 1000 << " sec / ";
        std::cout << t.elapsed() / 1000 / 60 << " min / ";
        std::cout << (t.elapsed() * 1000) / options.num_reads << " musec/read" << std::endl;
        std::cout << "num_mapped_reads " << options.num_mapped_reads << "/" << options.num_reads << " ("
                  << (options.num_mapped_reads * 100.0) / options.num_reads << "%)" << std::endl;
    }

    return 0;
}

template <typename FulgorIndex, typename QueryReader>
int pseudoalign_cobs_worker(FulgorIndex const& index, QueryReader& query_reader,
                            std::ofstream& out_file, std::mutex& out_mut,
                            const double threshold, ps_options& options) {
    cobs_output_buffer<FulgorIndex> output_buffer(index, out_file, out_mut);
    std::vector<pseudoalignment_match> matches;

    auto qg = query_reader.get_query_group();
    while (qg.refill()) {
        while (qg.has_next()) {
            typename QueryReader::query_t query;
            qg.value(query);

            index.pseudoalign_threshold_union(query.seq, matches, threshold);
            options.increment_processed_reads();
            if (!matches.empty()) options.increment_mapped_reads();
            output_buffer.write(query.name, matches);

            matches.clear();
            qg.next();
        }
    }

    return 0;
}

template <typename FulgorIndex, typename QueryReader>
int pseudoalign_cobs_orchestrator(FulgorIndex& index, QueryReader& query_reader,
                                  std::ofstream& out_file, const double threshold,
                                  ps_options& options) {
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> t;
    t.start();

    const uint64_t num_threads = options.num_threads;
    assert(num_threads >= 2);

    if (options.verbose) essentials::logger("*** START: pseudoalignment");
    std::vector<std::thread> workers;
    workers.reserve(num_threads);
    std::mutex out_mut;
    for (uint64_t i = 1; i != num_threads; ++i) {
        workers.push_back(std::thread([&index, &query_reader, &out_file, &out_mut, threshold,
                                       &options]() {
            pseudoalign_cobs_worker(index, query_reader, out_file, out_mut, threshold, options);
        }));
    }

    for (auto& w : workers) w.join();

    t.stop();
    if (options.verbose) essentials::logger("*** DONE: pseudoalignment");

    if (options.verbose) {
        std::cout << "processed " << options.num_reads << " reads" << std::endl;
        std::cout << "elapsed = " << t.elapsed() << " millisec / ";
        std::cout << t.elapsed() / 1000 << " sec / ";
        std::cout << t.elapsed() / 1000 / 60 << " min / ";
        std::cout << (t.elapsed() * 1000) / options.num_reads << " musec/read" << std::endl;
        std::cout << "num_mapped_reads " << options.num_mapped_reads << "/" << options.num_reads
                  << " (" << (options.num_mapped_reads * 100.0) / options.num_reads << "%)"
                  << std::endl;
    }

    return 0;
}

template <typename FulgorIndex, typename Formatter>
void fetch_and_deduplicate_sets(const std::string& query_filename,
                                Formatter& output_formatter,
                                std::string& tmp_filename,
                                FulgorIndex& index,
                                ps_options& options) {
    auto output_buffer = output_formatter.buffer();
    if (options.verbose) essentials::logger("*** START: fetching color set ids");

    std::ofstream tmp_file(tmp_filename, std::ios::binary);
    auto query_filenames = std::vector({query_filename});
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(query_filenames, options.num_threads,
                                                             options.num_threads - 1);
    rparser.start();
    std::vector<std::thread> workers;
    std::mutex outfile_mut, iomut;

    constexpr int32_t buff_thresh = 50;
    std::atomic<uint64_t> num_fetched_reads = 0;
    auto fetch = [&rparser, &index, &tmp_file, &outfile_mut, &iomut, &num_fetched_reads, &options] () {
        uint32_t buff_size = 0;
        std::vector<uint32_t> color_set_ids;
        std::stringstream ss;

        auto rg = rparser.getReadGroup();
        while (rparser.refill(rg)) {
            uint32_t read_id = rg.chunk_frag_offset().frag_idx;

            for (auto const& record: rg) {
                index.fetch_color_set_ids(record.seq, color_set_ids);

                buff_size += 1;

                ss.write(reinterpret_cast<char*>(&read_id), sizeof(read_id));
                uint32_t num_color_sets = static_cast<uint32_t>(color_set_ids.size());
                ss.write(reinterpret_cast<char*>(&num_color_sets), sizeof(num_color_sets));
                if (num_color_sets > 0) {
                    ss.write(reinterpret_cast<char*>(color_set_ids.data()), num_color_sets * sizeof(color_set_ids[0]));
                }

                color_set_ids.clear();
                if (options.verbose && num_fetched_reads > 0 && ++num_fetched_reads % 1000000 == 0) {
                    iomut.lock();
                    std::cout << "fetched " << num_fetched_reads << " reads" << std::endl;
                    iomut.unlock();
                }
                if (buff_size > buff_thresh) {
                    std::string outs = ss.str();
                    ss.str("");
                    outfile_mut.lock();
                    tmp_file.write(outs.data(), outs.size());
                    outfile_mut.unlock();
                    buff_size = 0;
                }
                ++read_id;
            }
        }
        if (buff_size > 0) {
            std::string outs = ss.str();
            ss.str("");
            outfile_mut.lock();
            tmp_file.write(outs.data(), outs.size());
            outfile_mut.unlock();
            buff_size = 0;
        }
    };

    for (uint64_t i = 1; i < options.num_threads; ++i) {
        workers.push_back(std::thread(fetch));
    }
    for (auto& w : workers) w.join();
    rparser.stop();
    tmp_file.close();

    if (options.verbose) essentials::logger("*** DONE: fetching color set ids");
    if (options.verbose) essentials::logger("*** START: deduplicating queries");

    std::ifstream ifile(tmp_filename, std::ios::binary);
    std::vector<std::vector<uint32_t>> queries;
    queries.reserve(num_fetched_reads);

    std::vector<uint32_t> tmp;
    uint32_t read_num = 0;
    while (ifile.read(reinterpret_cast<char*>(&read_num), sizeof(read_num))) {
        uint32_t num_colors = 0;
        ifile.read(reinterpret_cast<char*>(&num_colors), sizeof(num_colors));

        tmp.resize(num_colors + 1);
        tmp[0] = read_num;
        ifile.read(reinterpret_cast<char*>(&tmp[1]), num_colors * sizeof(num_colors));
        if (tmp.size() > 1) {
            queries.push_back(tmp);
            tmp.clear();
        } else {
            // just write out the unmapped reads here
            output_buffer.write(read_num, {});
        }
    }

    if (queries.empty()) { return; }

    std::sort(queries.begin(), queries.end(),
    [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) -> bool {
      return std::lexicographical_compare(a.begin() + 1, a.end(), b.begin() + 1,
                                          b.end());
    });

    auto curr = queries.begin();
    auto next = curr++;
    size_t identical_lists = 0;
    uint64_t identical_sizes = 0;
    while (next < queries.end()) {
        if ((curr->size() == next->size()) and
            std::equal(curr->begin() + 1, curr->end(), next->begin() + 1)) {
            identical_sizes += next->size() - 1;
            next->resize(1);  // retain only the id.
            ++identical_lists;
            ++next;
        } else {
            curr = next;
            ++next;
        }
    }

    std::cerr << "number of identical lists = " << identical_lists << " (skipping "
              << identical_sizes << " set ids)\n";
    ifile.close();
    std::ofstream ofile(tmp_filename, std::ios::trunc | std::ios::binary);
    for (auto& query : queries) {
        uint32_t s = query.size();
        ofile.write(reinterpret_cast<char*>(&s), sizeof(s));
        if (s > 0) {
            ofile.write(reinterpret_cast<char*>(query.data()), sizeof(query[0]) * s);
        }
    }
    ofile.close();

    if (options.verbose) essentials::logger("*** DONE: deduplicating queries");
}

int pseudoalign(int argc, char** argv) {
    std::vector<std::string> normalized_args;
    normalized_args.reserve(argc);
    for (int i = 0; i != argc; ++i) {
        normalized_args.emplace_back(argv[i]);
        if (normalized_args.back() == "--threshold") normalized_args.back() = "-r";
    }
    std::vector<char*> normalized_argv;
    normalized_argv.reserve(argc);
    for (auto& arg : normalized_args) normalized_argv.push_back(arg.data());

    cmd_line_parser::parser parser(argc, normalized_argv.data());

    parser.add("index_filename", "The Fulgor index filename.", "-i", true);
    parser.add("query_filename", "Query filename in FASTA/FASTQ format (optionally gzipped).", "-q",
               true);
    parser.add("output_filename",
               "File where output will be written. You can specify \"/dev/stdout\" to write "
               "output to stdout. In this case, it is also recommended to use the --verbose flag "
               "to avoid printing status messages to stdout.",
               "-o", true);
    parser.add("num_threads", "Number of threads (default is 1).", "-t", false);
    parser.add("verbose", "Verbose output during query (default is false).", "--verbose", false,
               true);
    parser.add("threshold",
               "Threshold for threshold_union algorithm. It must be a float in [0.0,1.0].", "-r",
               false);
    parser.add("cobs", "Write threshold-union output in the COBS-like format expected by Phylign.",
               "--cobs", false, true);
    parser.add("deduplicate", "Removes duplicate queries before pseudoalignment (default is false)."
               " Only works on Full-Intersection. Creates a temporary file in the executable's directory.",
               "--deduplicate", false, true);
    parser.add("format", "Format of the output file. Must either ascii, binary, compressed"
               " (default is ascii).", "--format", false);
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto query_filename = parser.get<std::string>("query_filename");
    auto output_filename = parser.get<std::string>("output_filename");

    bool deduplicate = parser.get<bool>("deduplicate");
    const bool cobs_output = parser.get<bool>("cobs");
    auto output_format = parser.parsed("format") ? parser.get<std::string>("format") : "ascii";

    uint64_t num_threads = 1;
    if (parser.parsed("num_threads")) num_threads = parser.get<uint64_t>("num_threads");
    if (num_threads == 1) {
        num_threads += 1;
        std::cerr
            << "1 thread was specified, but an additional thread will be allocated for parsing"
            << std::endl;
    }

    double threshold = constants::invalid_threshold;
    if (parser.parsed("threshold")) {
        threshold = parser.get<double>("threshold");
        if (threshold < 0.0 or threshold > 1.0) {
            std::cerr << "threshold must be a float in [0.0,1.0]" << std::endl;
            return 1;
        }
    }

    auto ps_alg = pseudoalignment_algorithm::FULL_INTERSECTION;
    if (threshold != constants::invalid_threshold) {
        if (deduplicate) {
            cerr << "Deduplication not available for threshold < 1.0. Remove --deduplicate flag." << std::endl;
            return 1;
        }
        ps_alg = pseudoalignment_algorithm::THRESHOLD_UNION;
    }
    if (cobs_output && ps_alg != pseudoalignment_algorithm::THRESHOLD_UNION) {
        std::cerr << "--cobs requires --threshold" << std::endl;
        return 1;
    }
    if (cobs_output && output_format != "ascii") {
        std::cerr << "--cobs requires --format ascii" << std::endl;
        return 1;
    }

    bool verbose = parser.get<bool>("verbose");
    if (verbose) util::print_cmd(argc, argv);

    std::variant<hfur_index_t, mdfur_index_t, mfur_index_t, dfur_index_t> index;
    if (is_meta_diff(index_filename)) {
        index = mdfur_index_t();
    } else if (is_meta(index_filename)) {
        index = mfur_index_t();
    } else if (is_diff(index_filename)) {
        index = dfur_index_t();
    } else if (is_hybrid(index_filename)) {
        index = hfur_index_t();
    } else {
        std::cerr << "Wrong index filename supplied." << std::endl;
        return 1;
    }

    std::variant<std::monostate, psa_ascii_formatter, psa_binary_formatter, psa_compressed_formatter> formatter;
    if (output_format == "ascii") {
        formatter.emplace<psa_ascii_formatter>(output_filename);
    } else if (output_format == "binary") {
        formatter.emplace<psa_binary_formatter>(output_filename);
    } else if (output_format == "compressed") {
        formatter.emplace<psa_compressed_formatter>(output_filename);
    } else {
        std::cout << "Unknown output format. Supported formats: ascii, binary, compressed." << std::endl;
        return 1;
    }

    std::string tmp_filename = "queries.tmp";
    ps_options options(ps_alg, verbose, num_threads);

    if (verbose) {
        std::cout << "\n---------------------------------" << std::endl;
        std::cout << "[Index]     " << index_filename << std::endl;
        std::cout << "[Queries]   " << query_filename << std::endl;
        std::cout << "[Output]    " << output_filename << std::endl;
        std::cout << "[Algorithm] " << to_string(ps_alg, threshold) << (deduplicate ? "(dedup.)" : "")
                  << (cobs_output ? " [cobs]" : "") << std::endl;
        std::cout << "---------------------------------\n" << std::endl;
    }

    std::visit([&index_filename, &query_filename, &output_filename, &tmp_filename,
                deduplicate, num_threads, threshold, verbose, cobs_output, &options]
                      (auto&& index, auto&& formatter) {
        if (verbose) essentials::logger("*** START: loading the index");
        essentials::load(index, index_filename.c_str());
        if (verbose) essentials::logger("*** DONE: loading the index");

        if (verbose) essentials::logger("performing queries from file '" + query_filename + "'...");

        if constexpr (std::is_same_v<std::decay_t<decltype(formatter)>, psa_compressed_formatter>) {
            formatter.set_num_colors(index.num_colors());
        }
        if constexpr (!std::is_same_v<std::decay_t<decltype(formatter)>, std::monostate>) {
            if (cobs_output) {
                std::ofstream out(output_filename);
                if (!out) {
                    std::cerr << "could not open output file " << output_filename << std::endl;
                    return;
                }
                named_fastq_query_reader query_reader(query_filename, num_threads, index);
                pseudoalign_cobs_orchestrator(index, query_reader, out, threshold, options);
            } else if (deduplicate) {
                fetch_and_deduplicate_sets(query_filename, formatter, tmp_filename, index, options);
                preprocessed_query_reader query_reader(tmp_filename, num_threads);
                pseudoalign_orchestrator(index, query_reader, formatter, threshold, options);
            } else {
                fastq_query_reader query_reader(query_filename, num_threads, index);
                pseudoalign_orchestrator(index, query_reader, formatter, threshold, options);
            }
        }

    }, index, formatter);

    std::remove(tmp_filename.c_str());

    return 0;
}
