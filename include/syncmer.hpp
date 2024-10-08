#pragma once

#include "util.hpp"
#include "enumerator.hpp"

namespace minimizers {

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
                pair.priority = m_params.p.open_syncmer;
            } else if (tmer_p == 0 or tmer_p == w0) {
                pair.priority = m_params.p.closed_syncmer;
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

            // compare kmers by kmer hash
            m_enum_kmers.eat_with_priority(kmer, m_params.k, preference);

            // compare kmers by tmer hash
            // m_enum_kmers.eat_with_priority(kmer + tmer_p, m_params.t, preference);
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

}  // namespace minimizers