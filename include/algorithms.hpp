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

/*
    Each algorithm returns a position p in [0..w-1], corresponding
    to the position of the kmer selected as the window's fingerprint.
    Note: in case of ties, we return the *leftmost* kmer.
*/

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

template <typename Hasher>
struct syncmer {
    syncmer(parameters const& params)
        : m_params(params)
        , m_offset((params.k - params.t) / 2)
        , m_enum_tmers(params.k - params.t + 1, params.t, params.seed)
        , m_enum_kmers(params.w, params.k, params.seed) {}

    uint64_t sample(char const* window) const {
        const uint64_t w0 = m_params.k - m_params.t;
        enumerator<Hasher> enum_tmers(w0 + 1, m_params.t, m_params.seed);
        pair_t<typename Hasher::hash_type> min_pair{uint64_t(-1), typename Hasher::hash_type(-1)};
        uint64_t p = -1;

        for (uint64_t i = 0; i != m_params.w; ++i) {
            char const* kmer = window + i;
            bool clear = i == 0;  // first kmer
            enum_tmers.eat(kmer, clear);
            uint64_t tmer_p = enum_tmers.next();
            assert(tmer_p <= w0);

            auto hash = Hasher::hash(kmer, m_params.k, m_params.seed);
            pair_t<typename Hasher::hash_type> pair{m_params.p.kmer, hash};

            if (tmer_p % m_params.w == m_offset) {
                pair.preference = m_params.p.open_syncmer;
            } else if (tmer_p == 0 or tmer_p == w0) {
                pair.preference = m_params.p.closed_syncmer;
            }

            if (pair < min_pair) {
                min_pair = pair;
                p = i;
            }
        }

        assert(p < m_params.w);
        return p;
    }

    uint64_t sample(char const* window, bool clear) {
        const uint64_t w0 = m_params.k - m_params.t;
        for (uint64_t i = clear ? 0 : m_params.w - 1; i != m_params.w; ++i) {
            char const* kmer = window + i;
            m_enum_tmers.eat(kmer, i == 0);
            uint64_t tmer_p = m_enum_tmers.next();
            assert(tmer_p <= w0);

            uint64_t preference = m_params.p.kmer;
            if (tmer_p % m_params.w == m_offset) {
                preference = m_params.p.open_syncmer;
            } else if (tmer_p == 0 or tmer_p == w0) {
                preference = m_params.p.closed_syncmer;
            }

            m_enum_kmers.eat_with_preference(kmer, m_params.k, preference);
        }

        uint64_t p = m_enum_kmers.next();
        assert(p < m_params.w);
        return p;
    }

private:
    parameters m_params;
    uint64_t m_offset;
    priority m_priority;
    enumerator<Hasher> m_enum_tmers;
    enumerator<Hasher> m_enum_kmers;
};

template <typename Anchor>
struct mod_sampling {
    mod_sampling(const uint64_t w, const uint64_t k, Anchor const& a)
        : m_w(w), m_k(k), m_M_w(fastmod::computeM_u32(m_w)), m_a(a) {}

    uint64_t sample(char const* window) const {
        uint64_t p = m_a.sample(window);
        return fastmod::fastmod_u32(static_cast<uint32_t>(p), m_M_w, m_w);  // p % w
    }

    uint64_t sample(char const* window, bool clear) {
        uint64_t p = m_a.sample(window, clear);
        return fastmod::fastmod_u32(static_cast<uint32_t>(p), m_M_w, m_w);  // p % w
    }

    uint64_t w() const { return m_w; }
    uint64_t k() const { return m_k; }

private:
    uint64_t m_w, m_k, m_M_w;
    Anchor m_a;
};

}  // namespace minimizers
