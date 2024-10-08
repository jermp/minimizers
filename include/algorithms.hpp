#pragma once

#include <cmath>
#include <numbers>
#include <complex>
#include <iomanip>

#include "external/fastmod/fastmod.h"

#include "util.hpp"
#include "enumerator.hpp"

namespace minimizers {

double redundancy_in_density_in_perc(const double density, const double lower_bound) {
    return (density / lower_bound - 1) * 100.0;
}

double redundancy_in_density_as_factor(const double density, const double lower_bound) {
    return density / lower_bound;
}

bool is_not_forward(const uint64_t k, const uint64_t w, const uint64_t t) {
    assert(w >= 2);
    assert(t <= k);
    /*
        We know a scheme is *not* foward when there exist x and y
        such that x mod w + 1 < y mod w, where x and y are the
        positions of the smallest t-mer in window i and i-1 respectively,
        for some i > 0. So we derive: x mod w < w - 2.
        Since x is in [0..l - t] = [0..w + k - 1 - t], then x is at most
        w + k - t - 1, i.e., w + k - t - 1 mod w < w - 2.

        All possible backward jumps (y mod w, x mod w), of length y-x-1,
        are for y in [x+1..w-1].

        Note: in math, we would write (k - t - 1) mod w < w - 2,
        but here we always sum w to avoid having to take the
        modulo of negative integers when t = k.
    */
    return (w + k - t - 1) % w < w - 2;
}

/* This ignores (asymptotic) lower order terms. */
double closed_form_density(std::string const& scheme_name,                        //
                           const uint64_t k, const uint64_t w, const uint64_t t)  //
{
    if (scheme_name == "miniception") {
        return 1.67 / w;  // at most
    } else if (scheme_name == "mod_sampling") {
        bool ok = (w + k - 1 - t) % w == w - 1;
        double correction = ok ? 0 : floor(1.0 + double(k - 1.0 - t) / w) / (w + k - t);
        return double(floor(1.0 + double(k - t - 1.0) / w) + 2.0 - correction) / (w + k - t + 1.0);
    } else {
        throw std::runtime_error("unknown scheme name");
    }
}

/*
    Each algorithm returns a position p in [0..w-1], corresponding
    to the position of the kmer selected as the window's fingerprint.
    Note: in case of ties, we return the *leftmost* kmer.
*/

template <typename Hasher>
struct mod_sampling {
    static std::string name() { return "mod_sampling"; }

    mod_sampling(uint64_t w, uint64_t k, uint64_t t, uint64_t seed)
        : m_w(w), m_k(k), m_t(t), m_seed(seed), m_enum_tmers(w + k - t, t, seed) {
        m_M_w = fastmod::computeM_u32(m_w);
    }

    uint64_t sample(char const* window) const {
        const uint64_t num_tmers = (m_w + m_k - 1) - m_t + 1;
        uint64_t p = -1;
        typename Hasher::hash_type min_hash(-1);
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
        //     for (; i < pos; ++i) { std::cout << "="; }
        //     std::cout << "|";
        //     if (p < pos + m_k) {
        //         for (; i < p; ++i) { std::cout << "="; }
        //         for (; i < p + m_t; ++i) { std::cout << "*"; }
        //         for (; i < pos + m_k; ++i) { std::cout << "="; }
        //         std::cout << "|";
        //     } else {
        //         for (; i < pos + m_k; ++i) { std::cout << "="; }
        //         std::cout << "|";
        //         for (; i < p; ++i) { std::cout << "="; }
        //         for (; i < p + m_t; ++i) { std::cout << "*"; }
        //     }
        //     for (; i < m_w + m_k - 1; ++i) { std::cout << "="; }
        //     std::cout << std::endl;
        // }

        return pos;
    }

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
struct closed_syncmer {
    static std::string name() { return "closed-syncmer"; }

    closed_syncmer(uint64_t w, uint64_t k, uint64_t t, uint64_t seed)
        : m_w(w)
        , m_k(k)
        , m_t(t)
        , m_seed(seed)
        , m_enum_tmers(k - t + 1, t, seed)
        , m_enum_kmers(w, k, seed) {}

    uint64_t sample(char const* window) const {
        const uint64_t w0 = m_k - m_t;
        enumerator<Hasher> enum_tmers(w0 + 1, m_t, m_seed);
        pair_t<typename Hasher::hash_type> min_pair{2, typename Hasher::hash_type(-1)};
        uint64_t p = -1;
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            bool clear = i == 0;  // first kmer
            enum_tmers.eat(kmer, clear);
            uint64_t tmer_p = enum_tmers.next();
            assert(tmer_p <= w0);
            auto hash = Hasher::hash(kmer, m_w, m_k, m_seed);
            pair_t<typename Hasher::hash_type> pair{1, hash};
            /* Prefer closed syncmers first, then kmers. */
            if (tmer_p == 0 or tmer_p == w0) pair.preference = 0;
            if (pair < min_pair) {
                min_pair = pair;
                p = i;
            }
        }
        assert(p < m_w);
        return p;
    }

    uint64_t sample(char const* window, bool clear) {
        const uint64_t w0 = m_k - m_t;
        for (uint64_t i = clear ? 0 : m_w - 1; i != m_w; ++i) {
            char const* kmer = window + i;
            m_enum_tmers.eat(kmer, i == 0);
            uint64_t tmer_p = m_enum_tmers.next();
            assert(tmer_p <= w0);
            uint64_t preference = 1;
            if (tmer_p == 0 or tmer_p == w0) preference = 0;
            m_enum_kmers.eat_with_preference(kmer, m_k, preference);
            // m_enum_kmers.eat_with_preference(kmer + tmer_p, m_t, preference);
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
struct open_syncmer {
    static std::string name() { return "open-syncmer"; }

    open_syncmer(uint64_t w, uint64_t k, uint64_t t, uint64_t seed)
        : m_w(w)
        , m_k(k)
        , m_t(t)
        , m_seed(seed)

        /* one could try other positions but (k-t)/2 is optimal for conservation */
        , m_p((m_k - m_t) / 2)

        , m_enum_tmers(k - t + 1, t, seed)
        , m_enum_kmers(w, k, seed) {}

    uint64_t sample(char const* window) const {
        const uint64_t w0 = m_k - m_t;
        enumerator<Hasher> enum_tmers(w0 + 1, m_t, m_seed);
        pair_t<typename Hasher::hash_type> min_pair{2, typename Hasher::hash_type(-1)};
        uint64_t p = -1;
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            bool clear = i == 0;  // first kmer
            enum_tmers.eat(kmer, clear);
            uint64_t tmer_p = enum_tmers.next();
            assert(tmer_p <= w0);
            auto hash = Hasher::hash(kmer, m_w, m_k, m_seed);
            pair_t<typename Hasher::hash_type> pair{1, hash};
            /* Prefer open syncmers first, then kmers. */
            if (tmer_p == m_p) pair.preference = 0;
            if (pair < min_pair) {
                min_pair = pair;
                p = i;
            }
        }
        assert(p < m_w);
        return p;
    }

    uint64_t sample(char const* window, bool clear) {
        const uint64_t w0 = m_k - m_t;
        for (uint64_t i = clear ? 0 : m_w - 1; i != m_w; ++i) {
            char const* kmer = window + i;
            m_enum_tmers.eat(kmer, i == 0);
            uint64_t tmer_p = m_enum_tmers.next();
            assert(tmer_p <= w0);
            uint64_t preference = 1;
            if (tmer_p == m_p) preference = 0;
            m_enum_kmers.eat_with_preference(kmer, m_k, preference);
            // m_enum_kmers.eat_with_preference(kmer + tmer_p, m_t, preference);
        }
        uint64_t p = m_enum_kmers.next();
        assert(p < m_w);
        return p;
    }

private:
    uint64_t m_w, m_k, m_t, m_seed;
    uint64_t m_p;  // if smallest tmer in kmer begins at position p, then it is an open syncmer
    enumerator<Hasher> m_enum_tmers;
    enumerator<Hasher> m_enum_kmers;
};

template <typename Hasher>
struct open_closed_syncmer {
    static std::string name() { return "open-closed-syncmer"; }

    open_closed_syncmer(uint64_t w, uint64_t k, uint64_t t, uint64_t seed)
        : m_w(w)
        , m_k(k)
        , m_t(t)
        , m_seed(seed)

        , m_p(((m_k - m_t) % m_w) / 2)

        , m_enum_tmers(k - t + 1, t, seed)
        , m_enum_kmers(w, k, seed) {}

    uint64_t sample(char const* window) const {
        const uint64_t w0 = m_k - m_t;
        enumerator<Hasher> enum_tmers(w0 + 1, m_t, m_seed);
        pair_t<typename Hasher::hash_type> min_pair{2, typename Hasher::hash_type(-1)};
        uint64_t p = -1;
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            bool clear = i == 0;  // first kmer
            enum_tmers.eat(kmer, clear);
            uint64_t tmer_p = enum_tmers.next();
            assert(tmer_p <= w0);

            // Warning: always hash kmers

            auto hash = Hasher::hash(kmer, m_w, m_k, m_seed);
            pair_t<typename Hasher::hash_type> pair{2, hash};

            if (tmer_p % m_w == m_p) {
                pair.preference = 0;
            } else if (tmer_p == 0 or tmer_p == w0) {
                pair.preference = 1;
            }

            if (pair < min_pair) {
                min_pair = pair;
                p = i;
            }
        }
        assert(p < m_w);
        return p;
    }

    uint64_t sample(char const* window, bool clear) {
        const uint64_t w0 = m_k - m_t;
        for (uint64_t i = clear ? 0 : m_w - 1; i != m_w; ++i) {
            char const* kmer = window + i;
            m_enum_tmers.eat(kmer, i == 0);
            uint64_t tmer_p = m_enum_tmers.next();
            assert(tmer_p <= w0);

            uint64_t preference = 2;
            if (tmer_p % m_w == m_p) {
                preference = 0;
            } else if (tmer_p == 0 or tmer_p == w0) {
                preference = 1;
            }

            // hash tmers rather than kmers when preference is 0
            m_enum_kmers.eat_with_preference(kmer + (preference == 0 ? tmer_p : 0),  //
                                             preference == 0 ? m_t : m_k,            //
                                             preference);

            // // always break ties by kmer hash
            // m_enum_kmers.eat_with_preference(kmer, m_k, preference);

            // // always break ties by tmer hash
            // m_enum_kmers.eat_with_preference(kmer + tmer_p, m_t, preference);
        }
        uint64_t p = m_enum_kmers.next();
        assert(p < m_w);
        return p;
    }

private:
    uint64_t m_w, m_k, m_t, m_seed;
    uint64_t m_p;  // if smallest tmer in kmer begins at position p, then it is an open syncmer
    enumerator<Hasher> m_enum_tmers;
    enumerator<Hasher> m_enum_kmers;
};

/*
    Miniception samples the smallest closed syncmer from a window.
*/
template <typename Hasher>
struct miniception {
    static std::string name() { return "miniception"; }

    miniception(uint64_t w, uint64_t k, uint64_t t, uint64_t seed) : m_alg(w, k, t, seed) {}

    uint64_t sample(char const* window) const { return m_alg.sample(window); }

    uint64_t sample(char const* window, bool clear) { return m_alg.sample(window, clear); }

private:
    closed_syncmer<Hasher> m_alg;
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

    uint64_t sample(char const* window) const {
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

uint8_t char_remap[256];

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
/// x        y>=x+sigma/2+1     y = leftmost max in row of x; must be created as some previous z+v,
/// v<=sigma-1. (z cannot be s0.)
///
/// p >= z+sigma/2+1>=y-(sigma-1)+sigma/2+1>=x+sigma/2+1-(sigma-1)+sigma/2+1=x+3 => must be created.
/// cannot be created by s0.

/// Return whether the kmer is in the UHS, and the random kmer order.
template <typename Hasher>
struct rotational_orig_hasher {
    using hash_type = uhs_hash<Hasher>;

    static hash_type hash(char const* kmer, const uint64_t w, const uint64_t k,
                          const uint64_t seed) {
        constexpr uint64_t sigma = 4;
        bool in_uhs = true;
        uint64_t sum0 = 0;
        for (uint64_t pos = 0; pos < k; pos += w) sum0 += char_remap[int(kmer[pos])];

        for (uint64_t j = 1; j != w; ++j) {
            uint64_t sumj = 0;
            for (uint64_t pos = j; pos < k; pos += w) sumj += char_remap[int(kmer[pos])];
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
        char_remap[int('A')] = 0;
        char_remap[int('C')] = 1;
        char_remap[int('T')] = 2;
        char_remap[int('G')] = 3;
    }

    uint64_t sample(char const* window) const {
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

// TODO: Global variables are ugly and ideally should be replaced by member variables on the Hasher
// objects.
vector<long double> sines;
vector<std::complex<long double>> roots;
long double pi = std::numbers::pi_v<long double>;

// The pseudocode from the original paper.
// We intentionally ignore the 0 case.
bool is_decycling_original(char const* kmer, const uint64_t k) {
    long double im = 0;
    for (uint64_t i = 0; i != k; ++i) { im += sines[i] * kmer[i]; }
    long double im_rot = 0;
    for (uint64_t i = 0; i != k; ++i) { im_rot += sines[i + 1] * kmer[i]; }
    // std::cerr << "Im: " << im << " im_rot: " << im_rot <<
    // std::endl;
    return im > 0 and im_rot <= 0;
}

// Same method, but using complex numbers.
// This is only different due to floating point rounding errors, e.g. when imag(x)=0.
// The original method has ever so slightly better density.
bool is_decycling_arg_pos(char const* kmer, const uint64_t k) {
    std::complex<long double> x = 0;
    for (uint64_t i = 0; i != k; ++i) { x += roots[i] * (long double)kmer[i]; }
    auto a = std::arg(x);
    // std::cerr << "arg:    " << a << "\nthresh: " << pi - 2 * pi / k << endl;
    return pi - 2 * pi / k < a;
}

// Use angle around 0 instead of around pi.
//
// This is the first negative instead of first positive rotation.
// That should be equivalent since it's basically using the D-tilde.
//
// FIXME: This is around 1% worse than the versions above. I do not understand why.
bool is_decycling_arg_neg(char const* kmer, const uint64_t k) {
    std::complex<long double> x = 0;
    for (uint64_t i = 0; i != k; ++i) { x += roots[i] * (long double)kmer[i]; }
    auto a = std::arg(x);
    return -2 * pi / k < a and a <= 0;
}

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
        sines.reserve(m_k + 1);
        for (uint64_t i = 0; i != m_k; ++i) { sines.push_back(std::sin(2 * pi * i / m_k)); }
        roots.reserve(m_k + 1);
        for (uint64_t i = 0; i != m_k; ++i) {
            roots.push_back(std::exp(std::complex<long double>(0, 2 * pi * i / m_k)));
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
        sines.reserve(m_k + 1);
        for (uint64_t i = 0; i != m_k; ++i) { sines.push_back(std::sin(2 * pi * i / m_k)); }
        roots.reserve(m_k + 1);
        for (uint64_t i = 0; i != m_k; ++i) {
            roots.push_back(std::exp(std::complex<long double>(0, 2 * pi * i / m_k)));
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
