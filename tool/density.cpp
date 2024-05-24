#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <sstream>

#include "external/cmd_line_parser/include/parser.hpp"
#include "external/gz/zip_stream.hpp"
#include "external/gz/zip_stream.cpp"

#include "hasher.hpp"
#include "algorithms.hpp"

template <typename T>
static inline void do_not_optimize_away(T&& value) {
    asm volatile("" : "+r"(value));
}

using namespace minimizers;

template <typename Alg>
void run(std::istream& is,                                                           //
         const uint64_t k, const uint64_t w, const uint64_t t, const uint64_t seed,  //
         const bool bench, const bool stream)                                        //
{
    typedef std::chrono::high_resolution_clock clock_type;
    const uint64_t l = w + k - 1;  // num. symbols in window

    std::string sequence;
    uint64_t pos = 0;
    uint64_t num_sampled_kmers = 0;
    uint64_t num_kmers = 0;
    uint64_t tot_num_windows = 0;

    bool is_forward = true;
    std::vector<uint64_t> positions;  // of sampled kmers

    Alg a(w, k, t, seed);

    auto duration = clock_type::duration::zero();

    while (util::appendline(is, sequence)) {
        if (sequence.size() == pos || sequence[pos] == '>' || sequence[pos] == ';') {
            // Empty line or new fasta entry.
            sequence.clear();
            continue;
        }
        if (sequence.size() < l) {
            // Not yet enough symbols for a window.
            // BUG: This shouldn't clear the sequence, and count kmers appropriately.
            sequence.clear();
            pos = 0;
            continue;
        }

        auto start = clock_type::now();

        const uint64_t num_windows = sequence.size() - l + 1;
        tot_num_windows += num_windows;
        for (uint64_t i = 0; i != num_windows; ++i) {
            char const* window = sequence.data() + i;
            uint64_t p = stream ? a.sample(window, i == 0) : a.sample(window);
            assert(p >= 0 and p < w);
            do_not_optimize_away(p);
            if (!bench) positions.push_back(i + p);  // absolute + relative
        }

        auto stop = clock_type::now();
        duration += stop - start;

        // Check forwardness outside of times loop.
        for (uint64_t i = 0; i + 1 < positions.size(); i++) {
            if (positions[i] > positions[i + 1]) {
                is_forward = false;
                break;
            }
        }

        // BUG: This is wrong when a sequence spans multiple lines.
        num_kmers += sequence.size() - k + 1;

        // Copy last l-1 symbols to the beginning of the buffer and drop the rest.
        std::copy(sequence.data() + sequence.size() - (l - 1), sequence.data() + sequence.size(),
                  sequence.data());
        sequence.resize(l - 1);
        pos = sequence.size();
    }

    if (bench) {
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
        std::cout << "tot. time: " << elapsed / 1000 << " [sec]" << std::endl;
        std::cout << "avg. per window: " << (elapsed * 1000000) / tot_num_windows << " [nanosec]"
                  << std::endl;
    } else {
        // BUG: when a sequence spans multiple lines, each line is counted again from 0.
        std::sort(positions.begin(), positions.end());
        num_sampled_kmers =
            std::distance(positions.begin(), std::unique(positions.begin(), positions.end()));

        std::stringbuf buffer;
        std::ostream os(&buffer);

        // std::cout << "  is_forward = " << (is_forward ? "YES" : "NO")
        //           << " (expected = " << (is_not_forward(k, w, t) ? "NO)" : "YES)") << std::endl;

        double density = double(num_sampled_kmers) / num_kmers;

        // os << "w,k,t,density,density_from_formula,is_forward\n";
        os << w << ',' << k << ',' << t << ',' << density << ',';

        std::cout << "  num_sampled_kmers = " << num_sampled_kmers << std::endl;
        std::cout << "  num_kmers = " << num_kmers << std::endl;
        std::cout << "  density = " << density << std::endl;
        std::cout << "  " << redundancy_in_density_as_factor(density, 1.0 / w)
                  << "X away from lower bound 1/w = " << 1.0 / w << std::endl;

        try {
            double density_from_formula = closed_form_density(Alg::name(), k, w, t);
            std::cout << "  calculation using closed formulas:\n";
            std::cout << "     density = " << density_from_formula << std::endl;
            std::cout << "     " << redundancy_in_density_as_factor(density_from_formula, 1.0 / w)
                      << "X away from lower bound 1/w = " << 1.0 / w << std::endl;
            os << density_from_formula << ',';
        } catch (std::runtime_error const& /*e*/) { os << ','; }

        os << (is_forward ? "YES" : "NO");

        std::cerr << buffer.str() << std::endl;
    }
}

template <typename Alg>
void run(std::string const& input_filename,  //
         const uint64_t k, const uint64_t w, const uint64_t t, const uint64_t seed,
         const bool bench, const bool stream)  //
{
    std::ifstream is(input_filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + input_filename + "'");
    if (util::ends_with(input_filename, ".gz")) {
        zip_istream zis(is);
        run<Alg>(zis, k, w, t, seed, bench, stream);
    } else {
        run<Alg>(is, k, w, t, seed, bench, stream);
    }
    is.close();
}

template <typename Hasher>
void run(std::string const& input_filename, std::string const& alg,  //
         const uint64_t k, const uint64_t w, const uint64_t seed,    //
         const bool bench, const bool stream)                        //
{
    const uint64_t r = 4;

    if (alg == "minimizer") {
        const uint64_t t = k;
        // TODO: The performance of the random minimizer should be measured using a dedicated
        // implementation that avoids the overhead of the modulo operation.
        run<mod_sampling<Hasher>>(input_filename, k, w, t, seed, bench, stream);
    } else if (alg == "lr-minimizer") {
        if (k >= w + r) {
            const uint64_t t = k - w;
            run<mod_sampling<Hasher>>(input_filename, k, w, t, seed, bench, stream);
        } else {
            std::cerr << "k must be at least w+r" << std::endl;
        }
    } else if (alg == "mod-minimizer") {
        const uint64_t t = r + ((k - r) % w);
        run<mod_sampling<Hasher>>(input_filename, k, w, t, seed, bench, stream);
    } else if (alg == "miniception") {
        // Theorem 7 of the miniception paper proves an upper bound for k-w+1,
        // but in practice they use k-w.
        const uint64_t t = w < k ? std::max(k - w, r) : r;
        run<miniception<Hasher>>(input_filename, k, w, t, seed, bench, stream);
    } else if (alg == "rot-minimizer-alt") {
        const uint64_t t = -1;  // not used
        run<rotational_alt<Hasher>>(input_filename, k, w, t, seed, bench, stream);
    } else if (alg == "rot-minimizer-orig") {
        if (k % w == 0) {
            const uint64_t t = -1;  // not used
            run<rotational_orig<Hasher>>(input_filename, k, w, t, seed, bench, stream);
        } else {
            std::cerr << "w must divide k" << std::endl;
        }
    } else if (alg == "mod-sampling") {
        for (uint64_t t = 1; t <= k; ++t) {
            std::cout << "### t = " << t << "\n";
            std::cout << "### num_tmers = " << (w + k - 1) - t + 1 << std::endl;
            run<mod_sampling<Hasher>>(input_filename, k, w, t, seed, bench, stream);
        }
    } else if (alg == "decycling") {
        const uint64_t t = -1;  // not used
        run<decycling<Hasher>>(input_filename, k, w, t, seed, bench, stream);
    } else if (alg == "double-decycling") {
        const uint64_t t = -1;  // not used
        run<double_decycling<Hasher>>(input_filename, k, w, t, seed, bench, stream);
    } else {
        std::cerr << "Error: '" << alg << "' does not correspond to any method" << std::endl;
    }
}

void run(std::string const& input_filename, std::string const& alg,                          //
         const uint64_t k, const uint64_t w, const uint64_t hash_size, const uint64_t seed,  //
         const bool bench, const bool stream)                                                //
{
    if (hash_size == 64) {
        run<hasher64_type>(input_filename, alg, k, w, seed, bench, stream);
    } else if (hash_size == 128) {
        run<hasher128_type>(input_filename, alg, k, w, seed, bench, stream);
    }
}

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("input_filename", "A fasta file (gzipped or not).", "-i", true);
    parser.add("k", "K-mer length.", "-k", true);
    parser.add("w", "Window size.", "-w", true);
    parser.add("alg",
               "Sampling algorithm to use. Options are: 'minimizer', 'lr-minimizer', "
               "'mod-minimizer', 'miniception', 'mod-sampling'.",
               "-a", true);
    parser.add("hash_size",
               "Number of bits for hash function (either 64 or 128; default is " +
                   std::to_string(constants::default_hash_size) + ").",
               "-b", false);
    parser.add(
        "seed",
        "Seed for hash function; default is " + std::to_string(constants::default_seed) + ".", "-s",
        false);
    parser.add("bench", "Benchmark sampling speed (do not compute density).", "--bench", false,
               true);
    parser.add("stream", "Sample from a stream (stateful computation rather than stateless).",
               "--stream", false, true);

    if (!parser.parse()) return 1;

    const auto input_filename = parser.get<std::string>("input_filename");
    const auto k = parser.get<uint64_t>("k");
    const auto w = parser.get<uint64_t>("w");
    const auto alg = parser.get<std::string>("alg");

    uint64_t hash_size = constants::default_hash_size;
    if (parser.parsed("hash_size")) {
        hash_size = parser.get<uint64_t>("hash_size");
        if (hash_size != 64 and hash_size != 128) {
            std::cerr << "Error: hash_size must be 64 or 128" << std::endl;
            return 1;
        }
    }

    uint64_t seed = constants::default_seed;
    if (parser.parsed("seed")) seed = parser.get<uint64_t>("seed");
    bool bench = parser.get<bool>("bench");
    bool stream = parser.get<bool>("stream");

    run(input_filename, alg, k, w, hash_size, seed, bench, stream);

    return 0;
}
