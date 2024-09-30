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

template <typename Hasher>
struct minimizer {
    static std::string name() { return "minimizer"; }

    minimizer(uint64_t w, uint64_t k, uint64_t /* t */, uint64_t seed)
        : m_w(w), m_k(k), m_seed(seed), m_enum_kmers(w, k, seed) {}

    uint64_t sample(char const* window) const {
        uint64_t p = -1;
        typename Hasher::hash_type min_hash(-1);
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            auto hash = Hasher::hash(kmer, m_k, m_seed);
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
    enumerator<Hasher> m_enum_kmers;
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
            auto hash = Hasher::hash(kmer, m_k, m_seed);
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
        , m_p((k - t) / 2)

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
            auto hash = Hasher::hash(kmer, m_k, m_seed);
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

        , m_p((k - t) / 2)

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

            auto hash = Hasher::hash(kmer, m_k, m_seed);
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

            // // hash tmers rather than kmers when preference is 0
            // m_enum_kmers.eat_with_preference(kmer + (preference == 0 ? tmer_p : 0),  //
            //                                  preference == 0 ? m_t : m_k,            //
            //                                  preference);

            // always break ties by kmer hash
            m_enum_kmers.eat_with_preference(kmer, m_k, preference);

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

template <typename Anchor>
struct mod_sampling {
    static std::string name() { return "mod_sampling"; }

    mod_sampling(const uint64_t w, const uint64_t k, Anchor const& a)
        : m_w(w), m_k(k), m_M_w(fastmod::computeM_u32(m_w)), m_a(a) {}

    uint64_t sample(char const* window) const {
        uint64_t p = m_a.sample(window);
        return fastmod::fastmod_u32(static_cast<uint32_t>(p), m_M_w, m_w);  // p % m_w
    }

    uint64_t sample(char const* window, bool clear) {
        uint64_t p = m_a.sample(window, clear);
        return fastmod::fastmod_u32(static_cast<uint32_t>(p), m_M_w, m_w);  // p % m_w
    }

    uint64_t w() const { return m_w; }
    uint64_t k() const { return m_k; }

private:
    uint64_t m_w, m_k, m_M_w;
    Anchor m_a;
};

}  // namespace minimizers
