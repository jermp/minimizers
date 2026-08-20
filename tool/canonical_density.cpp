#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "external/cmd_line_parser/include/parser.hpp"
#include "canonical_schemes.hpp"
#include "hasher.hpp"

/*
    Density of the canonical schemes of the paper

        D. Cologni and G. E. Pibiri, "Canonical Schemes and Minimizers".

    NOTE ON NOTATION. This tool uses the (k, m) notation of the paper,
    hence the parameters -k and -m: a window is a single k-mer, which
    contains w = k - m + 1 m-mers, and a scheme selects one m-mer out of
    every k-mer. This differs from the (w, k) notation assumed elsewhere
    in this repository (tool/density.cpp), where a window consists of w
    consecutive k-mers; the two conventions match under the substitution
    k -> m and w -> k - m + 1.

    For each scheme the tool reports the quantities of Table 3 of the
    paper: the density (fraction of distinct sampled positions, as a ratio
    to d = 2/(w+1) too), whether the scheme is forward, the number of
    sampled positions, and the number of super-k-mers (maximal runs of
    consecutive k-mers sampling the same position).
*/

using namespace minimizers;

template <typename Alg>
void run(std::string const& sequence, Alg alg) {
    const uint64_t k = alg.k();
    const uint64_t m = alg.m();
    const uint64_t w = alg.w();  // num. m-mers in a k-mer: w = k - m + 1

    if (sequence.length() < k) throw std::runtime_error("sequence is shorter than k");

    const uint64_t num_windows = sequence.length() - k + 1;  // num. k-mers
    const uint64_t num_mmers = sequence.length() - m + 1;    // num. sampleable positions

    std::vector<uint64_t> positions;  // of sampled m-mers, in window order
    positions.reserve(num_windows);
    for (uint64_t i = 0; i != num_windows; ++i) {
        char const* window = sequence.data() + i;
        uint64_t p = alg.sample(window);
        assert(p < w);
        positions.push_back(i + p);  // absolute = window start + relative
    }

    /* Forwardness: the sampled position never decreases as the window slides. */
    bool is_forward = true;
    for (uint64_t i = 0; i + 1 < positions.size(); ++i) {
        if (positions[i] > positions[i + 1]) {
            is_forward = false;
            break;
        }
    }

    /* Super-k-mers: maximal runs of consecutive windows sampling the same position. */
    uint64_t num_super_kmers = num_windows == 0 ? 0 : 1;
    for (uint64_t i = 1; i != positions.size(); ++i) {
        if (positions[i] != positions[i - 1]) num_super_kmers += 1;
    }

    std::sort(positions.begin(), positions.end());
    const uint64_t num_sampled_positions =
        std::distance(positions.begin(), std::unique(positions.begin(), positions.end()));

    const double d = 2.0 / (w + 1);  // density of the random minimizer, to leading order
    const double density = static_cast<double>(num_sampled_positions) / num_mmers;
    const double skm_rate = static_cast<double>(num_super_kmers) / num_windows;

    std::cout << "  k = " << k << ", m = " << m << ", w = k - m + 1 = " << w << std::endl;
    std::cout << "  is_forward = " << (is_forward ? "YES" : "NO") << std::endl;
    std::cout << "  num_sampled_positions = " << num_sampled_positions << std::endl;
    std::cout << "  num_super_kmers = " << num_super_kmers << std::endl;
    std::cout << "  num_positions = " << num_mmers << std::endl;
    std::cout << "  num_windows = " << num_windows << std::endl;
    std::cout << "  density = " << density << " (" << density / d
              << " x d, with d = 2/(w+1) = " << d << ")" << std::endl;
    std::cout << "  super_kmer_rate = " << skm_rate << " (" << skm_rate / d << " x d)" << std::endl;
}

template <typename Hasher>
void run(std::string const& input_filename, std::string const& alg_name,  //
         const uint64_t k, const uint64_t m, const uint64_t seed)         //
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

    if (alphabet_size % 2 != 0) {
        std::cerr << "Error: alphabet size must be even, so that complementation "
                     "(c -> alphabet_size - 1 - c) is fixed-point-free"
                  << std::endl;
        return;
    }

    canonical_parameters params(k, m, alphabet_size, seed);

    if (alg_name == "LvR") {
        run(sequence, lvr<Hasher>(params));
    } else if (alg_name == "RhL") {
        run(sequence, rhl<Hasher>(params));
    } else if (alg_name == "RR") {
        run(sequence, rr<Hasher>(params));
    } else if (alg_name == "LL") {
        run(sequence, ll<Hasher>(params));
    } else if (alg_name == "CK") {
        run(sequence, ck<Hasher>(params));
    } else {
        std::cerr << "Error: '" << alg_name << "' does not correspond to any method" << std::endl;
    }
}

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("input_filename",
               "A binary file encoding a sequence. It must be created with the program "
               "'generate_random_sequence' (with an even alphabet size).",
               "-i", true);
    parser.add("k",
               "K-mer length: the window. Note that this tool uses the (k, m) notation of the "
               "paper 'Canonical Schemes and Minimizers', not the (w, k) notation assumed "
               "elsewhere in this repository: a window is a k-mer, containing w = k - m + 1 "
               "m-mers.",
               "-k", true);
    parser.add("m", "Minimizer (m-mer) length, with m < k.", "-m", true);
    parser.add("alg",
               "Canonical sampling scheme to use. Options are: 'LvR', 'RhL', 'RR', 'LL', 'CK'.",
               "-a", true);
    parser.add("hash_size",
               "Number of bits for hash function (either 64 or 128; default is " +
                   std::to_string(constants::default_hash_size) + ").",
               "-b", false);
    parser.add(
        "seed",
        "Seed for hash function; default is " + std::to_string(constants::default_seed) + ".", "-s",
        false);

    if (!parser.parse()) return 1;

    const auto input_filename = parser.get<std::string>("input_filename");
    const auto k = parser.get<uint64_t>("k");
    const auto m = parser.get<uint64_t>("m");
    const auto alg = parser.get<std::string>("alg");

    if (m >= k) {
        std::cerr << "Error: it must be m < k" << std::endl;
        return 1;
    }

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

    if (hash_size == 64) {
        run<hasher64_type>(input_filename, alg, k, m, seed);
    } else {
        run<hasher128_type>(input_filename, alg, k, m, seed);
    }

    return 0;
}
