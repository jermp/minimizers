#pragma once

#include <cmath>

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
        return 1.67 / w;
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
        m_M_w = fastmod::computeM_u64(m_w);
    }

    uint64_t sample(char const* window) const {
        const uint64_t num_tmers = (m_w + m_k - 1) - m_t + 1;
        uint64_t p = 0;
        typename Hasher::hash_type min_hash(-1);
        for (uint64_t i = 0; i != num_tmers; ++i) {
            char const* tmer = window + i;
            auto hash = Hasher::hash(tmer, m_t, m_seed);
            if (hash < min_hash) {
                min_hash = hash;
                p = i;
            }
        }
        return fastmod::fastmod_u64(p, m_M_w, m_w);  // p % m_w
    }

    uint64_t sample(char const* window, bool clear) {
        m_enum_tmers.eat(window, clear);
        uint64_t p = m_enum_tmers.next();
        return fastmod::fastmod_u64(p, m_M_w, m_w);  // p % m_w
    }

private:
    uint64_t m_w, m_k, m_t, m_seed;
    __uint128_t m_M_w;
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
        uint64_t p = 0;
        typename Hasher::hash_type min_hash(-1);
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            bool clear = i == 0;  // first kmer
            enum_tmers.eat(kmer, clear);
            uint64_t tmer_p = enum_tmers.next();
            assert(tmer_p >= 0 and tmer_p <= w0);
            if (tmer_p == 0 or tmer_p == w0) {  // context is charged
                auto hash = Hasher::hash(kmer, m_k, m_seed);
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

struct rotational {
    static std::string name() { return "rotational"; }

    rotational(uint64_t w, uint64_t k, uint64_t /*t*/, uint64_t /*seed*/) : m_w(w), m_k(k) {}

    uint64_t sample(char const* window) {
        uint64_t p = 0;
        uint64_t max = 0;
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            uint64_t sum = 0;
            for (uint64_t j = 0; j < m_k; j += m_w) sum += kmer[j];
            if (sum > max) {
                max = sum;
                p = i;
            }
        }
        assert(p < m_w);
        return p;
    }

    uint64_t sample(char const* window, bool /*clear*/) {
        // Warning: not implemented...
        return sample(window);
    }

private:
    uint64_t m_w, m_k;
};

}  // namespace minimizers