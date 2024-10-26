#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <sstream>
#include <numeric>  // for std::accumulate

#include "external/cmd_line_parser/include/parser.hpp"
#include "hasher.hpp"
#include "mod_sampling.hpp"
#include "syncmer.hpp"
#include "double_decycling.hpp"

template <typename T>
static inline void do_not_optimize_away(T&& value) {
    asm volatile("" : "+r"(value));
}

using namespace minimizers;

template <typename Alg>
double run(std::string const& sequence, Alg alg, const bool bench, const bool stream) {
    typedef std::chrono::high_resolution_clock clock_type;

    uint64_t num_sampled_kmers = 0;
    bool is_forward = true;
    std::vector<uint64_t> positions;  // of sampled kmers

    auto duration = clock_type::duration::zero();
    auto start = clock_type::now();

    const uint64_t w = alg.w();
    const uint64_t k = alg.k();
    const uint64_t l = w + k - 1;  // num. symbols in window
    const uint64_t num_windows = sequence.length() - l + 1;
    const uint64_t num_kmers = sequence.length() - k + 1;

    for (uint64_t i = 0; i != num_windows; ++i) {
        char const* window = sequence.data() + i;
        uint64_t p = stream ? alg.sample(window, i == 0) : alg.sample(window);
        assert(p < w);
        do_not_optimize_away(p);
        if (!bench) {
            uint64_t absolute_pos = i + p;  // absolute + relative
            positions.push_back(absolute_pos);
        }
    }

    auto stop = clock_type::now();
    duration += stop - start;

    /* Check forwardness outside main loop. */
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

    os << w << ',' << k << ',' << density << ',';

    std::cout << "  num_sampled_kmers = " << num_sampled_kmers << std::endl;
    std::cout << "  num_kmers = " << num_kmers << std::endl;
    std::cout << "  num_windows = " << num_windows << std::endl;
    std::cout << "  density = " << density << std::endl;
    std::cout << "  " << util::redundancy_in_density_as_factor(density, 1.0 / w)
              << "X away from lower bound 1/w = " << 1.0 / w << std::endl;

    os << (is_forward ? "YES" : "NO");
    std::cerr << buffer.str() << std::endl;

    return density;
}

template <typename Hasher>
void run(std::string const& input_filename, std::string const& alg_name,  //
         const uint64_t k, const uint64_t w, const uint64_t seed,         //
         const bool bench, const bool stream)                             //
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

    typedef syncmer<Hasher> anchor_type;

    /*
        Choose r based on alphabet size.
        (There are just examples based on our experiments:
         we do not aim at being exhaustive here.)
    */
    uint64_t r = alphabet_size <= 4 ? 4 : 1;
    uint64_t s = k > w ? std::max(k - w, r) : r;
    uint64_t t = k;

    if (util::begins_with(alg_name, "mod")) {
        t = r + ((k - r) % w);
        s = r;
    }

    priority p;
    if (alg_name == "DD" or alg_name == "mod-DD") {
        p = {0, 0, 0};
        parameters params(w + k - t, t, s, seed, p);
        typedef double_decycling<Hasher> anchor_type;
        anchor_type anchor(params);
        mod_sampling<anchor_type> alg(w, k, anchor);
        run(sequence, alg, bench, stream);
        return;
    } else if (alg_name == "M" or alg_name == "mod-M") {
        p = {0, 0, 0};
    } else if (alg_name == "C" or alg_name == "mod-C") {
        p = {1, 0, 2};
    } else if (alg_name == "O" or alg_name == "mod-O") {
        p = {1, 2, 0};
    } else if (alg_name == "OC" or alg_name == "mod-OC") {
        p = {2, 1, 0};
    } else {
        std::cerr << "Error: '" << alg_name << "' does not correspond to any method" << std::endl;
        return;
    }

    parameters params(w + k - t, t, s, seed, p);
    anchor_type anchor(params);
    mod_sampling<anchor_type> alg(w, k, anchor);
    run(sequence, alg, bench, stream);
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
               "Sampling algorithm to use. Options are: 'M', 'C', 'O', 'OC', 'DD', "
               "'mod-M', 'mod-C', 'mod-O', 'mod-OC', 'mod-DD'.",
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