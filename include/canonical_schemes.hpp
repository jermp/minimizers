#pragma once

#include <cassert>
#include <cstdint>
#include <cstring>
#include <vector>

#include "util.hpp"

namespace minimizers {

/*
    Canonical sampling schemes from the paper

        D. Cologni and G. E. Pibiri, "Canonical Schemes and Minimizers".

    NOTE ON NOTATION. This file, and the tool tool/canonical_density.cpp,
    use the (k, m) notation of the paper: a window is a single k-mer, which
    contains w = k - m + 1 m-mers, and a scheme selects one m-mer out of
    every k-mer. The rest of this repository uses the (w, k) notation
    instead, where a window consists of w consecutive k-mers. The two
    conventions match under the substitution k -> m and w -> k - m + 1:
    a window there is a k-mer here.

    Writing rc(s) for the reverse complement of a string s,
    canonical(s) := min_lex{s, rc(s)} for its canonical form,
    x_i for the m-mer of the k-mer x at position i in [0, k - m], and
    MINI(x) for the ordinary (leftmost) random minimizer of x under a
    hash function h, the five schemes (Definitions 7 and 8 of the paper) are:

        LvR(x) := min_h   { canonical(x_0), ..., canonical(x_{k-m}) }
        RhL(x) := min_lex { MINI(x), MINI(rc(x)) }
        RR(x)  := min_h   { min_h{x_0, rc(x_0)}, ..., min_h{x_{k-m}, rc(x_{k-m})} }
        LL(x)  := min_lex { canonical(x_0), ..., canonical(x_{k-m}) }
        CK(x)  := MINI(canonical(x))

    All five are canonical: x and rc(x) receive the same signature
    (Proposition 5 and Definition 7 of the paper).

    LvR, RR and LL rank each position by a key that depends only on the
    reverse-complement class of the m-mer sitting there ("schemes on C"),
    so a tie among positions carrying the same class must be broken in a
    mirror-equivariant way: no leftmost rule can be, because the mirror
    i -> (k - m) - i reverses the order of the tied positions
    (Proposition 23). These three schemes therefore break ties with the
    centre-based tie-break of Proposition 24: among the tied positions,
    take the one closest to the centre (k - m)/2 of the window, and when
    two are equally close (they are then mirror images i and (k - m) - i)
    take the smaller index if x <= rc(x) and the larger otherwise. This
    rule is forward for all k and m, and mirror-equivariant for odd k.
    The same tie-break is implemented in SSHash
    (https://github.com/jermp/sshash, branch alt-canon-mode,
    file include/util.hpp).

    The alphabet is the one of tool/generate_random_sequence.cpp: symbols
    are the raw bytes 0, 1, ..., alphabet_size - 1, compared bytewise, and
    the complement of symbol c is alphabet_size - 1 - c, a fixed-point-free
    involution for even alphabet sizes. For alphabet_size = 4 this is DNA
    under the encoding A=0, C=1, G=2, T=3, with A <-> T and C <-> G.
*/

struct canonical_parameters {
    canonical_parameters() : canonical_parameters(0, 0, 0, 0) {}
    canonical_parameters(uint64_t k, uint64_t m, uint64_t alphabet_size, uint64_t seed)
        : k(k), m(m), alphabet_size(alphabet_size), seed(seed) {}
    uint64_t k, m, alphabet_size, seed;
};

template <typename Hasher>
struct canonical_scheme {
    typedef typename Hasher::hash_type hash_type;

    canonical_scheme(canonical_parameters const& params) : m_params(params), m_rc(params.k) {
        assert(params.m > 0 and params.m < params.k);
        assert(params.alphabet_size % 2 == 0);  // complementation must be fixed-point-free
    }

    uint64_t k() const { return m_params.k; }
    uint64_t m() const { return m_params.m; }
    uint64_t w() const { return m_params.k - m_params.m + 1; }  // num. m-mers in a k-mer

protected:
    canonical_parameters m_params;
    std::vector<char> m_rc;  // reverse complement of the current k-mer

    void fill_rc(char const* kmer) {
        const uint64_t k = m_params.k;
        const char sigma_minus_1 = static_cast<char>(m_params.alphabet_size - 1);
        for (uint64_t i = 0; i != k; ++i) m_rc[i] = sigma_minus_1 - kmer[k - 1 - i];
    }

    /* rc(x_i) = rc(x)_{(k-m)-i}, by the mirror identity (Eq. (1) of the paper). */
    char const* rc_mmer(const uint64_t i) const {
        assert(i <= m_params.k - m_params.m);
        return m_rc.data() + (m_params.k - m_params.m - i);
    }

    /* The canonical form of the m-mer of x at position i. */
    char const* canonical_mmer(char const* kmer, const uint64_t i) const {
        char const* fwd = kmer + i;
        char const* rc = rc_mmer(i);
        return std::memcmp(fwd, rc, m_params.m) <= 0 ? fwd : rc;
    }

    /* Negative if x < rc(x), zero if x = rc(x), positive if x > rc(x). */
    int compare_strands(char const* kmer) const {
        return std::memcmp(kmer, m_rc.data(), m_params.k);
    }

    hash_type hash_mmer(char const* mmer) const {
        return Hasher::hash(mmer, m_params.m, m_params.seed);
    }

    /*
        The ordinary (leftmost) random minimizer of the w m-mers of the
        string `strand`: positions are ranked by the pair (hash, m-mer),
        so that two distinct m-mers ranked equally by the hash are
        compared lexicographically (as required in Section 2 of the paper),
        and the leftmost position carrying the winning m-mer is returned.
    */
    void strand_minimizer(char const* strand, uint64_t& pos, char const*& mmer) const {
        pos = 0;
        mmer = strand;
        hash_type min_hash = hash_mmer(strand);
        for (uint64_t i = 1; i != w(); ++i) {
            char const* candidate = strand + i;
            hash_type hash = hash_mmer(candidate);
            if (hash < min_hash or
                (hash == min_hash and std::memcmp(candidate, mmer, m_params.m) < 0)) {
                min_hash = hash;
                mmer = candidate;
                pos = i;
            }
        }
    }

    /*
        The minimizer of a scheme on C: positions are ranked by the pair
        (key_hash(i), canonical m-mer at i), and ties among the positions
        carrying the winning class are broken with the centre-based rule
        of Proposition 24. `key_hash` must depend on position i only
        through the class of the m-mer sitting there.
    */
    template <typename KeyHash>
    uint64_t sample_on_C(char const* kmer, KeyHash key_hash) const {
        const uint64_t m = m_params.m;

        /* Find the smallest key (hash, canonical m-mer). */
        char const* min_mmer = canonical_mmer(kmer, 0);
        hash_type min_hash = key_hash(0);
        bool tie = false;
        for (uint64_t i = 1; i != w(); ++i) {
            char const* mmer = canonical_mmer(kmer, i);
            hash_type hash = key_hash(i);
            if (hash < min_hash) {
                min_hash = hash;
                min_mmer = mmer;
                tie = false;
            } else if (hash == min_hash) {
                int cmp = std::memcmp(mmer, min_mmer, m);
                if (cmp < 0) {
                    min_mmer = mmer;
                    tie = false;
                } else if (cmp == 0) {
                    tie = true;
                }
            }
        }

        /*
            Tie-break of Proposition 24: among the tied positions, all
            carrying the same class, take the one closest to the centre
            (k-m)/2 of the window; the doubled distance |2i - (k-m)|
            avoids the half-integer centre of an even window. When two
            positions are equally close they are mirror images (lo, hi)
            with lo + hi = k - m: take lo if x <= rc(x), hi otherwise.
        */
        uint64_t lo = 0, hi = 0;
        const uint64_t two_c = m_params.k - m;
        uint64_t best_dist = uint64_t(-1);
        for (uint64_t i = 0; i != w(); ++i) {
            if (key_hash(i) == min_hash and
                std::memcmp(canonical_mmer(kmer, i), min_mmer, m) == 0) {
                const uint64_t dist = 2 * i > two_c ? 2 * i - two_c : two_c - 2 * i;
                if (dist < best_dist) {
                    best_dist = dist;
                    lo = i;
                    hi = i;
                } else if (dist == best_dist) {
                    hi = i;
                }
                if (!tie) break;  // single occurrence of the winning class: done
            }
        }
        assert(lo == hi or lo + hi == two_c);
        return (lo == hi or compare_strands(kmer) <= 0) ? lo : hi;
    }
};

/*
    LvR: replace each m-mer by its canonical form and take the ordinary
    random minimizer of the resulting sequence. Forward, mirror-equivariant
    (for odd k), of density d = 2/(w+1): Theorem 15 of the paper.
*/
template <typename Hasher>
struct lvr : canonical_scheme<Hasher> {
    lvr(canonical_parameters const& params) : canonical_scheme<Hasher>(params) {}

    uint64_t sample(char const* kmer) {
        this->fill_rc(kmer);
        return this->sample_on_C(
            kmer, [&](const uint64_t i) { return this->hash_mmer(this->canonical_mmer(kmer, i)); });
    }
};

/*
    RR: rank each position by min{h(x_i), h(rc(x_i))}, a key symmetric in
    the pair {x_i, rc(x_i)}. Forward, mirror-equivariant (for odd k), of
    density d = 2/(w+1): Proposition 21 of the paper.
*/
template <typename Hasher>
struct rr : canonical_scheme<Hasher> {
    typedef typename Hasher::hash_type hash_type;

    rr(canonical_parameters const& params) : canonical_scheme<Hasher>(params) {}

    uint64_t sample(char const* kmer) {
        this->fill_rc(kmer);
        return this->sample_on_C(kmer, [&](const uint64_t i) {
            hash_type fwd = this->hash_mmer(kmer + i);
            hash_type rc = this->hash_mmer(this->rc_mmer(i));
            return fwd < rc ? fwd : rc;
        });
    }
};

/*
    LL: take the lexicographically smallest canonical m-mer. Forward and
    mirror-equivariant (for odd k); its density has no closed form, being
    a property of the string rather than a constant (Section 5 of the
    paper). Kraken and the KMC line instantiate this shape, with scrambled
    or restricted orders in place of plain lexicographic.
*/
template <typename Hasher>
struct ll : canonical_scheme<Hasher> {
    ll(canonical_parameters const& params) : canonical_scheme<Hasher>(params) {}

    uint64_t sample(char const* kmer) {
        this->fill_rc(kmer);
        /* A constant key: positions are compared by canonical m-mer alone. */
        return this->sample_on_C(kmer, [](const uint64_t) { return uint64_t(0); });
    }
};

/*
    RhL: minimize each strand under h, then keep the lexicographically
    smaller of the two winners. This is the canonical scheme of SSHash.
    Not forward; it opens super-k-mers at rate (4/3)d to leading order,
    because it selects under one order and compares under another:
    Theorem 19 of the paper.
*/
template <typename Hasher>
struct rhl : canonical_scheme<Hasher> {
    rhl(canonical_parameters const& params) : canonical_scheme<Hasher>(params) {}

    uint64_t sample(char const* kmer) {
        this->fill_rc(kmer);
        uint64_t fwd_pos, rc_pos;
        char const *fwd_min, *rc_min;
        this->strand_minimizer(kmer, fwd_pos, fwd_min);
        this->strand_minimizer(this->m_rc.data(), rc_pos, rc_min);
        int cmp = std::memcmp(fwd_min, rc_min, this->m());
        /*
            When MINI(x) = MINI(rc(x)) the comparison does not determine
            which occurrence of the common winner is sampled and
            mirror-equivariance fails (footnote 2 of the paper); we keep
            the forward strand's position in that case.
        */
        if (cmp <= 0) return fwd_pos;
        /* The winner sits at position rc_pos of rc(x), i.e. at the
           mirror position (k-m) - rc_pos of x. */
        return (this->k() - this->m()) - rc_pos;
    }
};

/*
    CK: canonicalize the whole k-mer, then take its ordinary minimizer,
    sampled at the position of x at which it sits. The trivial baseline:
    canonical and mirror-equivariant for every tie-break, but not forward,
    with density -> 2d and Theta(1) super-k-mers per context: Theorem 20
    of the paper.
*/
template <typename Hasher>
struct ck : canonical_scheme<Hasher> {
    ck(canonical_parameters const& params) : canonical_scheme<Hasher>(params) {}

    uint64_t sample(char const* kmer) {
        this->fill_rc(kmer);
        int cmp = this->compare_strands(kmer);
        char const* xstar = cmp <= 0 ? kmer : this->m_rc.data();
        uint64_t pos;
        char const* mmer;
        this->strand_minimizer(xstar, pos, mmer);
        /* If the canonical k-mer is rc(x), position pos of rc(x) is
           position (k-m) - pos of x. */
        return cmp <= 0 ? pos : (this->k() - this->m()) - pos;
    }
};

}  // namespace minimizers
