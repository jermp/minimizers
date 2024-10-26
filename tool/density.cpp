#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <sstream>
#include <numeric>  // for std::accumulate

#include "external/cmd_line_parser/include/parser.hpp"
#include "hasher.hpp"
#include "algorithms.hpp"

template <typename T>
static inline void do_not_optimize_away(T&& value) {
    asm volatile("" : "+r"(value));
}

using namespace minimizers;

template <typename Alg>
double run(std::string const& sequence,                                                //
           const uint64_t k, const uint64_t w, const uint64_t t, const uint64_t seed,  //
           const bool bench, const bool stream)                                        //
{
    typedef std::chrono::high_resolution_clock clock_type;

    uint64_t num_sampled_kmers = 0;
    bool is_forward = true;
    std::vector<uint64_t> positions;  // of sampled kmers

    Alg a(w, k, t, seed);

    auto duration = clock_type::duration::zero();
    auto start = clock_type::now();

    const uint64_t l = w + k - 1;  // num. symbols in window
    const uint64_t num_windows = sequence.length() - l + 1;
    const uint64_t num_kmers = sequence.length() - k + 1;

    for (uint64_t i = 0; i != num_windows; ++i) {
        char const* window = sequence.data() + i;
        uint64_t p = stream ? a.sample(window, i == 0) : a.sample(window);
        assert(p < w);
        do_not_optimize_away(p);
        if (!bench) {
            uint64_t absolute_pos = i + p;  // absolute + relative
            positions.push_back(absolute_pos);
        }
    }

    auto stop = clock_type::now();
    duration += stop - start;

    // Check forwardness outside main loop.
    for (uint64_t i = 0; i + 1 < positions.size(); i++) {
        if (positions[i] > positions[i + 1]) {
            is_forward = false;
            break;
        }
    }

    if (bench) {
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
        std::cout << "tot. time: " << elapsed / 1000 << " [sec]" << std::endl;
        std::cout << "avg. per window: " << (elapsed * 1000000) / num_windows << " [nanosec]"
                  << std::endl;
        return 0.0;
    }
    std::sort(positions.begin(), positions.end());
    num_sampled_kmers =
        std::distance(positions.begin(), std::unique(positions.begin(), positions.end()));

    std::stringbuf buffer;
    std::ostream os(&buffer);

    std::cout << "  is_forward = " << (is_forward ? "YES" : "NO") << std::endl;

    double density = static_cast<double>(num_sampled_kmers) / num_kmers;

    os << w << ',' << k << ',' << t << ',' << density << ',';

    std::cout << "  num_sampled_kmers = " << num_sampled_kmers << std::endl;
    std::cout << "  num_kmers = " << num_kmers << std::endl;
    std::cout << "  num_windows = " << num_windows << std::endl;

    std::cout << "  density = " << density << std::endl;
    std::cout << "  " << redundancy_in_density_as_factor(density, 1.0 / w)
              << "X away from lower bound 1/w = " << 1.0 / w << std::endl;

    try {
        double density_from_formula = closed_form_density(Alg::name(), k, w, t);
        std::cout << "  calculation using closed formulas:\n";
        std::cout << "     density = " << density_from_formula << std::endl;
        std::cout << "     " << redundancy_in_density_as_factor(density_from_formula, 1.0 / w)
                  << "X away from lower bound 1/w = " << 1.0 / w << std::endl;
        // os << density_from_formula << ',';
    } catch (std::runtime_error const& /*e*/) { os << ','; }

    os << (is_forward ? "YES" : "NO");
    std::cerr << buffer.str() << std::endl;

    return density;
}

template <typename Hasher>
void run(std::string const& input_filename, std::string const& alg,  //
         const uint64_t k, const uint64_t w, const uint64_t seed,    //
         const bool bench, const bool stream)                        //
{
    std::ifstream is(input_filename.c_str(), std::ifstream::binary);
    if (!is.good()) throw std::runtime_error("error in opening the file '" + input_filename + "'");

    /* read input sequence */
    std::string sequence;
    uint64_t sequence_length = 0;
    uint64_t alphabet_size = 0;
    is.read(reinterpret_cast<char*>(&sequence_length), sizeof(sequence_length));
    is.read(reinterpret_cast<char*>(&alphabet_size), sizeof(alphabet_size));
    assert(alphabet_size <= 256);
    sequence.resize(sequence_length);
    is.read(sequence.data(), sequence_length);
    is.close();

    /*
        Choose r based on alphabet size.
        (There are just examples based on our experiments:
         we do not aim at being exhaustive here.)
    */
    uint64_t r = alphabet_size <= 4 ? 4 : 1;

    if (alg == "minimizer") {
        const uint64_t t = k;
        run<mod_sampling<Hasher>>(sequence, k, w, t, seed, bench, stream);
    }

    else if (alg == "lr-minimizer") {
        if (k >= w + r) {
            const uint64_t t = k - w;
            run<mod_sampling<Hasher>>(sequence, k, w, t, seed, bench, stream);
        } else {
            std::cerr << "k must be at least w+r" << std::endl;
        }
    }

    else if (alg == "closed-syncmer") {
        const uint64_t t = k > w ? std::max(k - w, r) : r;
        run<closed_syncmer<Hasher>>(sequence, k, w, t, seed, bench, stream);
    }

    else if (alg == "open-syncmer") {
        const uint64_t t = k > w ? std::max(k - w, r) : r;
        run<open_syncmer<Hasher>>(sequence, k, w, t, seed, bench, stream);
    }

    else if (alg == "open-closed-syncmer") {
        uint64_t t = r;
        run<open_closed_syncmer<Hasher>>(sequence, k, w, t, seed, bench, stream);
    }

    else if (alg == "mod-minimizer") {
        uint64_t t = r + ((k - r) % w);
        run<mod_sampling<Hasher>>(sequence, k, w, t, seed, bench, stream);
    }

    else if (alg == "miniception") {
        // Theorem 7 of the miniception paper proves an upper bound for k-w+1,
        // but in practice they use k-w.
        const uint64_t t = k > w ? std::max(k - w, r) : r;
        run<miniception<Hasher>>(sequence, k, w, t, seed, bench, stream);
    }

    else if (alg == "rot-minimizer-alt") {
        const uint64_t t = -1;  // not used
        run<rotational_alt<Hasher>>(sequence, k, w, t, seed, bench, stream);
    }

    else if (alg == "rot-minimizer-orig") {
        if (k % w == 0) {
            const uint64_t t = -1;  // not used
            run<rotational_orig<Hasher>>(sequence, k, w, t, seed, bench, stream);
        } else {
            std::cerr << "w must divide k" << std::endl;
        }
    }

    else if (alg == "mod-sampling") {
        for (uint64_t t = 1; t <= k; ++t) {
            std::cout << "### t = " << t << "\n";
            std::cout << "### num_tmers = " << (w + k - 1) - t + 1 << std::endl;
            run<mod_sampling<Hasher>>(sequence, k, w, t, seed, bench, stream);
        }
    }

    else if (alg == "decycling") {
        const uint64_t t = -1;  // not used
        run<decycling<Hasher>>(sequence, k, w, t, seed, bench, stream);
    }

    else if (alg == "double-decycling") {
        const uint64_t t = -1;  // not used
        run<double_decycling<Hasher>>(sequence, k, w, t, seed, bench, stream);
    }

    else {
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
    parser.add("input_filename",
               "A binary file encoding a sequence. It must be created with the program "
               "'generate_random_sequence'.",
               "-i", true);
    parser.add("k", "K-mer length.", "-k", true);
    parser.add("w", "Window size.", "-w", true);
    parser.add("alg",
               "Sampling algorithm to use. Options are: 'minimizer', 'lr-minimizer', "
               "'mod-minimizer', 'miniception', 'mod-sampling', 'rot-minimizer-orig', "
               "'rot-minimizer-alt', 'decycling', 'double-decycling', 'closed-syncmer', "
               "'open-syncmer', 'open-closed-syncmer'.",
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
