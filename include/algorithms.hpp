#pragma once

#include <cmath>
#include <numbers>
#include <complex>
#include <cassert>
#include <iostream>
#include <array>

#include "external/fastmod/fastmod.h"

#include "enumerator.hpp"

namespace minimizers {

namespace detail {

struct globals_t {
    std::vector<long double> sines;
    std::vector<std::complex<long double>> roots;
    long double pi = std::numbers::pi_v<long double>;
};

globals_t& get_globals();

}  // namespace detail

double redundancy_in_density_in_perc(double density, double lower_bound);
double redundancy_in_density_as_factor(double density, double lower_bound);
bool is_not_forward(const uint64_t k, const uint64_t w, const uint64_t t);

/*
    Each algorithm returns a position p in [0..w-1], corresponding
    to the position of the kmer selected as the window's fingerprint.
    Note: in case of ties, we return the *leftmost* kmer.
*/

double closed_form_density(std::string const& scheme_name, const uint64_t k, const uint64_t w,
                           const uint64_t t);

template <typename Hasher>
struct mod_sampling {
    static std::string name() { return "mod_sampling"; }

    mod_sampling(uint64_t w, uint64_t k, uint64_t t, uint64_t seed)
        : m_w(w), m_k(k), m_t(t), m_seed(seed), m_enum_tmers(w + k - t, t, seed) {
        m_M_w = fastmod::computeM_u32(m_w);
    }

    /// Sample from a single window.
    uint64_t sample(char const* window) const {
        const uint64_t num_tmers = (m_w + m_k - 1) - m_t + 1;
        uint64_t p = -1;
        typename Hasher::hash_type min_hash(-1);
        // Find the leftmost tmer with minimal hash.
        for (uint64_t i = 0; i != num_tmers; ++i) {
            char const* tmer = window + i;
            auto hash = Hasher::hash(tmer, m_w, m_t, m_seed);
            if (hash < min_hash) {
                min_hash = hash;
                p = i;
            }
        }
        assert(p < num_tmers);
        uint64_t pos = fastmod::fastmod_u32(p, m_M_w, m_w);  // p % m_w

        // if (p == pos) {
        //     uint64_t i = 0;
        //     for (; i != p; ++i) { std::cout << "="; }
        //     std::cout << "|";
        //     for (; i != p + m_t; ++i) { std::cout << "*"; }
        //     for (; i != p + m_k; ++i) { std::cout << "="; }
        //     std::cout << "|";
        //     for (; i != m_w + m_k - 1; ++i) { std::cout << "="; }
        //     std::cout << std::endl;
        // } else {
        //     assert(pos < p);
        //     uint64_t i = 0;
        //     for (; i != pos; ++i) { std::cout << "="; }
        //     std::cout << "|";
        //     for (; i != p; ++i) { std::cout << "="; }
        //     for (; i != p + m_t; ++i) { std::cout << "*"; }
        //     for (; i != pos + m_k; ++i) { std::cout << "="; }
        //     std::cout << "|";
        //     for (; i != m_w + m_k - 1; ++i) { std::cout << "="; }
        //     std::cout << std::endl;
        // }

        return pos;
    }

    /// Sample from a stream.
    /// If `clear`, this is the first call.
    uint64_t sample(char const* window, bool clear) {
        m_enum_tmers.eat(window, clear);
        uint64_t p = m_enum_tmers.next();
        return fastmod::fastmod_u32(p, m_M_w, m_w);  // p % m_w
    }

private:
    uint64_t m_w, m_k, m_t, m_seed;
    uint64_t m_M_w;
    enumerator<Hasher> m_enum_tmers;
};

template <typename Hasher>
struct miniception {
    static std::string name() { return "miniception"; }

    miniception(uint64_t w, uint64_t k, uint64_t t, uint64_t seed)
        : m_w(w)
        , m_k(k)
        , m_t(t)
        , m_seed(seed)
        , m_enum_tmers(k - t + 1, t, seed)
        , m_enum_kmers(w, k, seed) {}

    uint64_t sample(char const* window) const {
        const uint64_t w0 = m_k - m_t;
        enumerator<Hasher> enum_tmers(w0 + 1, m_t, m_seed);
        uint64_t p = -1;
        typename Hasher::hash_type min_hash(-1);
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            bool clear = i == 0;  // first kmer
            enum_tmers.eat(kmer, clear);
            uint64_t tmer_p = enum_tmers.next();
            assert(tmer_p >= 0 and tmer_p <= w0);
            if (tmer_p == 0 or tmer_p == w0) {  // context is charged
                auto hash = Hasher::hash(kmer, m_w, m_k, m_seed);
                if (hash < min_hash) {
                    min_hash = hash;
                    p = i;
                }
            }
        }
        assert(p < m_w);
        return p;
    }

    uint64_t sample(char const* window, bool clear) {
        for (uint64_t i = clear ? 0 : m_w - 1; i != m_w; ++i) {
            char const* kmer = window + i;
            m_enum_tmers.eat(kmer, i == 0);
            uint64_t tmer_p = m_enum_tmers.next();
            assert(tmer_p >= 0 and tmer_p <= m_k - m_t);
            if (tmer_p == 0 or tmer_p == m_k - m_t) {  // context is charged
                m_enum_kmers.eat(kmer);
            } else {
                m_enum_kmers.skip();
            }
        }
        uint64_t p = m_enum_kmers.next();
        assert(p < m_w);
        return p;
    }

private:
    uint64_t m_w, m_k, m_t, m_seed;
    enumerator<Hasher> m_enum_tmers;
    enumerator<Hasher> m_enum_kmers;
};

template <typename Hasher>
using rotational_alt_hash = std::pair<int64_t, typename Hasher::hash_type>;

/// Return the negative of the sum of characters in positions 0 mod w, so that
/// the kmer with max sum compares smallest.
template <typename Hasher>
struct rotational_alt_hasher {
    using hash_type = rotational_alt_hash<Hasher>;

    // TODO: This can be implemented in O(1) time by storing prefix sums and using a rolling hash.
    static hash_type hash(char const* kmer, const uint64_t w, const uint64_t k,
                          const uint64_t seed) {
        int64_t sum = 0;
        for (uint64_t j = 0; j < k; j += w) sum += kmer[j];
        return {-sum, Hasher::hash(kmer, w, k, seed)};
    }
};

/// Our own simpler and much faster version.
/// Sample the leftmost kmer with the largest sum of characters in positions 0 mod w.
/// This is equivalent to a mod_sampling with the rotational_hasher function.
template <typename Hasher>
struct rotational_alt {
    static std::string name() { return "rotational_alt"; }

    rotational_alt(uint64_t w, uint64_t k, uint64_t /*t*/, uint64_t seed)
        : m_w(w), m_k(k), m_seed(seed), m_enum_kmers(w, k, seed) {}

    uint64_t sample(char const* window) {
        uint64_t p = -1;
        rotational_alt_hash<Hasher> min_hash{-1, -1};
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            auto hash = rotational_alt_hasher<Hasher>::hash(kmer, m_w, m_k, m_seed);
            if (hash < min_hash) {
                min_hash = hash;
                p = i;
            }
        }
        assert(p < m_w);
        return p;
    }

    uint64_t sample(char const* window, bool clear) {
        m_enum_kmers.eat(window, clear);
        return m_enum_kmers.next();
    }

private:
    uint64_t m_w, m_k, m_seed;
    enumerator<rotational_alt_hasher<Hasher>> m_enum_kmers;
};

template <typename Hasher>
using uhs_hash = std::pair<uint8_t, typename Hasher::hash_type>;

/// Proof that for all j: sumj <= sum0 + sigma/2 is sufficient to guarantee the
/// existence of a sum0.
///
/// We have a grid
/// sum0    sumj
/// s0   s1 ... sw-1               s0+sigma/2+1
/// s1   s2 ... s0+x               s1+sigma/2+1
/// ...
/// sw-1 s0+x s1+y ... sw-1+z      s[w-1]+sigma/w+1
///
/// Suppose that for each row, there is a j such that sumj <= su0+sigma/2 does not hold.
/// I.e. there is a value at least ^
///
/// s0       x>=s0+sigma/2+1    x = leftmost max in the 1st row
///
/// .                        p=q+u
/// ..
/// z        p>=z+sigma/2+1          (implied by y below)
/// .                        y=z+v
/// ..
/// x        y>=x+sigma/2+1     y = leftmost max in row of x; must be created as some previous
/// z+v, v<=sigma-1. (z cannot be s0.)
///
/// p >= z+sigma/2+1>=y-(sigma-1)+sigma/2+1>=x+sigma/2+1-(sigma-1)+sigma/2+1=x+3 => must be
/// created. cannot be created by s0.

/// Return whether the kmer is in the UHS, and the random kmer order.
template <typename Hasher>
struct rotational_orig_hasher {
    using hash_type = uhs_hash<Hasher>;

    static hash_type hash(char const* kmer, const uint64_t w, const uint64_t k,
                          const uint64_t seed) {
        static constexpr uint64_t sigma = 4;
        static constexpr auto char_remap = []() -> std::array<char, 256> {
            std::array<char, 256> dst;

            dst[int('A')] = 0;
            dst[int('C')] = 1;
            dst[int('T')] = 2;
            dst[int('G')] = 3;

            return dst;
        }();

        bool in_uhs = true;
        uint64_t sum0 = 0;
        for (uint64_t pos = 0; pos < k; pos += w) sum0 += char_remap.at(int(kmer[pos]));

        for (uint64_t j = 1; j != w; ++j) {
            uint64_t sumj = 0;
            for (uint64_t pos = j; pos < k; pos += w) sumj += char_remap.at(int(kmer[pos]));
            // Assume alphabet size 4.
            // Instead of <=+sigma, we do <=+sigma-1,
            // since the max difference between two characters is actually
            // sigma-1, not sigma.
            // In fact, I have a sketch of a proof (see above) that sigma/2 is sufficient.
            // (We make sure to do the correct rounding in odd cases.)
            // if (!(sumj <= sum0 + (sigma + 1) / 2)) {
            if (!(sumj <= sum0 + sigma)) {
                in_uhs = false;
                break;
            }
        }

        return {in_uhs ? 0 : 1, Hasher::hash(kmer, w, k, seed)};
    }
};

/* Version faithful to the original description by Marcais et al. */
template <typename Hasher>
struct rotational_orig {
    static std::string name() { return "rotational_orig"; }

    rotational_orig(uint64_t w, uint64_t k, uint64_t /*t*/, uint64_t seed)
        : m_w(w), m_k(k), m_seed(seed), m_enum_kmers(w, k, seed) {
        assert(m_k % m_w == 0);
    }

    uint64_t sample(char const* window) {
        uint64_t p = -1;
        uhs_hash<Hasher> min_hash{-1, -1};
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            auto hash = rotational_orig_hasher<Hasher>::hash(kmer, m_w, m_k, m_seed);
            if (hash < min_hash) {
                min_hash = hash;
                p = i;
            }
        }
        if (min_hash.first != 0) {
            std::cerr << "Not a single kmer is in UHS!" << std::endl;
            std::exit(1);
        }
        assert(p < m_w);
        return p;
    }

    uint64_t sample(char const* window, bool clear) {
        m_enum_kmers.eat(window, clear);
        return m_enum_kmers.next();
    }

private:
    uint64_t m_w, m_k, m_seed;
    enumerator<rotational_orig_hasher<Hasher>> m_enum_kmers;
};

// The pseudocode from the original paper.
// We intentionally ignore the 0 case.
bool is_decycling_original(char const* kmer, const uint64_t k);

// Same method, but using complex numbers.
// This is only different due to floating point rounding errors, e.g. when imag(x)=0.
// The original method has ever so slightly better density.
bool is_decycling_arg_pos(char const* kmer, const uint64_t k);

// Use angle around 0 instead of around pi.
//
// This is the first negative instead of first positive rotation.
// That should be equivalent since it's basically using the D-tilde.
//
// FIXME: This is around 1% worse than the versions above. I do not understand why.
bool is_decycling_arg_neg(char const* kmer, const uint64_t k);

template <typename Hasher>
struct decycling_hasher {
    using hash_type = uhs_hash<Hasher>;

    // TODO: This can be implemented in O(1) using rolling embedding.
    static hash_type hash(char const* kmer, const uint64_t w, const uint64_t k,
                          const uint64_t seed) {
        // for testing
        // auto is_decycling_1 = is_decycling_original(kmer, k);
        auto is_decycling_2 = is_decycling_arg_pos(kmer, k);
        // auto is_decycling_3 = is_decycling_arg_neg(kmer, k);
        auto is_decycling = is_decycling_2;
        return {is_decycling ? 0 : 1, Hasher::hash(kmer, w, k, seed)};
    }
};

template <typename Hasher>
struct decycling {
    static std::string name() { return "decycling"; }

    decycling(uint64_t w, uint64_t k, uint64_t /*t*/, uint64_t seed)
        : m_w(w), m_k(k), m_seed(seed), m_enum_kmers(w, k, seed) {
        detail::get_globals().sines.reserve(m_k + 1);
        for (uint64_t i = 0; i != m_k; ++i) {
            detail::get_globals().sines.push_back(std::sin(2 * detail::get_globals().pi * i / m_k));
        }
        detail::get_globals().roots.reserve(m_k + 1);
        for (uint64_t i = 0; i != m_k; ++i) {
            detail::get_globals().roots.push_back(
                std::exp(std::complex<long double>(0, 2 * detail::get_globals().pi * i / m_k)));
        }
    }

    uint64_t sample(char const* window) {
        using T = std::pair<uint8_t, typename Hasher::hash_type>;
        uint64_t p = -1;
        T min_hash{-1, -1};
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            auto hash = decycling_hasher<Hasher>::hash(kmer, m_w, m_k, m_seed);
            if (hash < min_hash) {
                min_hash = hash;
                p = i;
            }
        }
        assert(p < m_w);
        return p;
    }

    uint64_t sample(char const* window, bool clear) {
        m_enum_kmers.eat(window, clear);
        return m_enum_kmers.next();
    }

private:
    uint64_t m_w, m_k, m_seed;
    enumerator<decycling_hasher<Hasher>> m_enum_kmers;
};

template <typename Hasher>
struct double_decycling_hasher {
    using hash_type = uhs_hash<Hasher>;

    // TODO: This can be implemented in O(1) using rolling embedding.
    static hash_type hash(char const* kmer, const uint64_t w, const uint64_t k,
                          const uint64_t seed) {
        // FIXME: Using _original instead of _pos gives slightly better density?
        bool is_decycling_pos = is_decycling_arg_pos(kmer, k);
        bool is_decycling_neg = is_decycling_arg_neg(kmer, k);
        if (is_decycling_pos) { return {0, Hasher::hash(kmer, w, k, seed)}; }
        if (is_decycling_neg) { return {1, Hasher::hash(kmer, w, k, seed)}; }
        return {2, Hasher::hash(kmer, w, k, seed)};
    }
};

template <typename Hasher>
struct double_decycling {
    static std::string name() { return "double_decycling"; }

    double_decycling(uint64_t w, uint64_t k, uint64_t /*t*/, uint64_t seed)
        : m_w(w), m_k(k), m_seed(seed), m_enum_kmers(w, k, seed) {
        detail::get_globals().sines.reserve(m_k + 1);
        for (uint64_t i = 0; i != m_k; ++i) {
            detail::get_globals().sines.push_back(std::sin(2 * detail::get_globals().pi * i / m_k));
        }
        detail::get_globals().roots.reserve(m_k + 1);
        for (uint64_t i = 0; i != m_k; ++i) {
            detail::get_globals().roots.push_back(
                std::exp(std::complex<long double>(0, 2 * detail::get_globals().pi * i / m_k)));
        }
    }

    uint64_t sample(char const* window) {
        using T = std::pair<uint8_t, typename Hasher::hash_type>;
        uint64_t p = -1;
        T min_hash{-1, -1};
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            auto hash = double_decycling_hasher<Hasher>::hash(kmer, m_w, m_k, m_seed);
            if (hash < min_hash) {
                min_hash = hash;
                p = i;
            }
        }
        assert(p < m_w);
        return p;
    }

    uint64_t sample(char const* window, bool clear) {
        m_enum_kmers.eat(window, clear);
        return m_enum_kmers.next();
    }

private:
    uint64_t m_w, m_k, m_seed;
    enumerator<double_decycling_hasher<Hasher>> m_enum_kmers;
};

}  // namespace minimizers
