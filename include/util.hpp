#pragma once

namespace minimizers {

struct priority {
    priority() : priority(0, 0, 0) {}
    priority(uint64_t kmer, uint64_t closed_syncmer, uint64_t open_syncmer)
        : kmer(kmer), closed_syncmer(closed_syncmer), open_syncmer(open_syncmer) {}
    uint64_t kmer;
    uint64_t closed_syncmer;
    uint64_t open_syncmer;
};

struct parameters {
    parameters() : parameters(0, 0, 0, 0, priority()) {}
    parameters(uint64_t w, uint64_t k, uint64_t t, uint64_t seed, priority p)
        : w(w), k(k), t(t), seed(seed), p(p) {}
    uint64_t w, k, t, seed;
    priority p;
};

namespace constants {
constexpr uint64_t default_hash_size = 128;
constexpr uint64_t default_seed = 1234567890;
}  // namespace constants

namespace util {

double redundancy_in_density_in_perc(const double density, const double lower_bound) {
    return (density / lower_bound - 1) * 100.0;
}

double redundancy_in_density_as_factor(const double density, const double lower_bound) {
    return density / lower_bound;
}

static void print_cmd(int argc, char** argv) {
    for (int i = 0; i != argc; ++i) std::cout << argv[i] << ' ';
    std::cout << std::endl;
}

static bool begins_with(std::string const& str, std::string const& pattern) {
    if (pattern.size() > str.size()) return false;
    return std::equal(pattern.begin(), pattern.end(), str.begin());
}

static bool ends_with(std::string const& str, std::string const& pattern) {
    if (pattern.size() > str.size()) return false;
    return std::equal(pattern.begin(), pattern.end(), str.end() - pattern.size());
}

}  // namespace util

}  // namespace minimizers
