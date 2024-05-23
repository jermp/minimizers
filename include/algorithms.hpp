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

/// Return the negative of the sum of characters in positions 0 mod w, so that
/// the kmer with max sum compares smallest.
struct rotational_hasher {
    typedef int64_t hash_type;

    rotational_hasher(uint64_t w) : m_w(w) {}

    static int64_t hash(char const* kmer, const uint64_t w, const uint64_t k,
                        const uint64_t /*seed*/) {
        int64_t sum = 0;
        for (uint64_t j = 0; j < k; j += w) sum += kmer[j];
        return -sum;
    }

private:
    uint64_t m_w;
};

/// Our own simpler and much faster version.
/// Sample the leftmost kmer with the largest sum of characters in positions 0 mod w.
/// This is equivalent to a mod_sampling with the rotational_hasher function.
struct rotational_alt {
    static std::string name() { return "rotational_alt"; }

    rotational_alt(uint64_t w, uint64_t k, uint64_t /*t*/, uint64_t /*seed*/) : m_w(w), m_k(k) {}

    uint64_t sample(char const* window) {
        uint64_t p = -1;
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

    uint64_t sample(char const* window, bool clear) {
        m_enum_kmers.eat(window, clear);
        auto p = m_enum_kmers.next();
        // FIXME: Test and then remove.
        auto q = sample(window);
        assert(p == q);
        return p;
    }

private:
    uint64_t m_w, m_k;
    enumerator<rotational_hasher> m_enum_kmers;
};

/* Version faithful to the original description by Marcais et al. */
template <typename Hasher>
struct rotational_orig {
    static std::string name() { return "rotational_orig"; }

    rotational_orig(uint64_t w, uint64_t k, uint64_t /*t*/, uint64_t seed)
        : m_w(w), m_k(k), m_seed(seed) {
        assert(m_k % m_w == 0);
    }

    uint64_t sample(char const* window) {
        using T = std::pair<uint8_t, typename Hasher::hash_type>;
        static std::vector<uint64_t> sum(m_w);
        uint64_t p = -1;
        T min_hash{-1, -1};
        for (uint64_t i = 0; i != m_w; ++i) {
            char const* kmer = window + i;
            std::fill(sum.begin(), sum.end(), 0);
            for (uint64_t j = 0; j != m_k; ++j) sum[j % m_w] += kmer[j];
            bool in_uhs = true;
            for (uint64_t j = 1; j != m_w; ++j) {
                if (!(sum[j] <= sum[0] + 4)) {  // assume alphabet of size 4
                    in_uhs = false;
                    break;
                }
            }
            T hash{in_uhs ? 0 : 1, Hasher::hash(kmer, m_k, m_seed)};
            if (hash < min_hash) { min_hash = hash p = i; }
        }
        assert(p < m_w);
        return p;
    }

    uint64_t sample(char const* window, bool /*clear*/) {
        // Warning: not implemented...
        return sample(window);
    }

private:
    uint64_t m_w, m_k, m_seed;
};

/// Decycling code from original paper.
/// TODO: double decycling
template <typename Hasher>
struct decycling {
    static std::string name() { return "decycling-1"; }

    decycling(uint64_t w, uint64_t k, uint64_t /*t*/, uint64_t seed)
        : m_w(w), m_k(k), m_seed(seed), pi(std::numbers::pi_v<long double>) {
        m_sines.reserve(m_k + 1);
        for (uint64_t i = 0; i != m_k + 1; ++i) { m_sines.push_back(std::sin(2 * pi * i / m_k)); }
        m_roots.reserve(m_k + 1);
        for (uint64_t i = 0; i != m_k + 1; ++i) {
            m_roots.push_back(std::exp(std::complex<long double>(0, 2 * pi * i / m_k)));
        }
        std::cerr << std::setprecision(20);
    }

    // The pseudocode from the original paper.
    // We intentionally ignore the 0 case.
    bool is_decycling_original(char const* kmer) {
        long double im = 0;
        for (uint64_t i = 0; i != m_k; ++i) { im += m_sines[i] * kmer[i]; }
        long double im_rot = 0;
        for (uint64_t i = 0; i != m_k; ++i) { im_rot += m_sines[i + 1] * kmer[i]; }
        // std::cerr << "Im: " << im << " im_rot: " << im_rot <<
        // std::endl;
        return im > eps and im_rot <= eps;
    }

    // Same method, but using complex numbers.
    // The only time this is different is in rare cases where im_rot=0,
    // but the original method above computes im_rot as just above 0 due to rounding errors.
    bool is_decycling_arg_same(char const* kmer) {
        std::complex<long double> x = 0;
        for (uint64_t i = 0; i != m_k; ++i) { x += m_roots[i] * (long double)kmer[i]; }
        auto a = std::arg(x);
        // std::cerr << "arg:    " << a << "\nthresh: " << pi - 2 * pi / m_k << endl;
        return pi - 2 * pi / m_k + eps < a;
    }

    // Use angle around 0 instead of around pi.
    //
    // This is the first negative instead of first positive rotation.
    // That should be equivalent since it's basically using the D-tilde.
    //
    // FIXME: This is around 0.5% worse than the versions above. I do not understand why.
    bool is_decycling_arg_equiv(char const* kmer) {
        std::complex<long double> x = 0;
        for (uint64_t i = 0; i != m_k; ++i) { x += m_roots[i] * (long double)kmer[i]; }
        auto a = std::arg(x);
        return -2 * pi / m_k + eps < a and a <= eps;
    }

    // More natural version that does anti-clockwise instead of clockwise negative rotation.
    //
    // FIXME: This is around 5% worse than the versions above. I do not understand why.
    // I expected this to be equivalent, since the rotation direction shouldn't matter?
    bool is_decycling_arg_simple(char const* kmer) {
        std::complex<long double> x = 0;
        for (uint64_t i = 0; i != m_k; ++i) { x += m_roots[i] * (long double)kmer[i]; }
        auto a = std::arg(x);
        return -eps <= a and a < 2 * pi / m_k - eps;
    }

    uint64_t sample(char const* window) {
        using T = std::pair<uint8_t, typename Hasher::hash_type>;
        uint64_t p = -1;
        T min_hash{-1, -1};
        for (uint64_t i = 0; i != m_w; ++i) {
            auto is_decycling_1 = is_decycling_original(window + i);
            auto is_decycling_2 = is_decycling_arg_same(window + i);
            auto is_decycling_3 = is_decycling_arg_equiv(window + i);
            auto is_decycling_4 = is_decycling_arg_simple(window + i);
            auto is_decycling = is_decycling_2;
            T hash{is_decycling ? 0 : 1, Hasher::hash(window + i, m_w, m_k, m_seed)};
            if (hash < min_hash) {
                min_hash = hash;
                p = i;
            }
        }
        assert(p < m_w);
        return p;
    }

    uint64_t sample(char const* window, bool /*clear*/) {
        // TODO
        return sample(window);
    }

private:
    long double eps = 1e-14;
    long double pi;
    uint64_t m_w, m_k, m_seed;
    vector<long double> m_sines;
    vector<std::complex<long double>> m_roots;
};

}  // namespace minimizers
