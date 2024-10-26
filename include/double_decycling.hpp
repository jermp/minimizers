#pragma once

#include <cmath>
#include <numbers>
#include <complex>

#include "util.hpp"
#include "enumerator.hpp"

namespace minimizers {

vector<long double> sines;
vector<std::complex<long double>> roots;
long double pi = std::numbers::pi_v<long double>;

template <typename Hasher>
using uhs_hash = std::pair<uint8_t, typename Hasher::hash_type>;

template <typename Hasher>
struct double_decycling_hasher {
    using hash_type = uhs_hash<Hasher>;
    static hash_type hash(char const* kmer, const uint64_t k, const uint64_t seed) {
        bool is_decycling_pos = is_decycling_arg_pos(kmer, k);
        bool is_decycling_neg = is_decycling_arg_neg(kmer, k);
        if (is_decycling_pos) { return {0, Hasher::hash(kmer, k, seed)}; }
        if (is_decycling_neg) { return {1, Hasher::hash(kmer, k, seed)}; }
        return {2, Hasher::hash(kmer, k, seed)};
    }

private:
    static bool is_decycling_arg_pos(char const* kmer, const uint64_t k) {
        std::complex<long double> x = 0;
        for (uint64_t i = 0; i != k; ++i) { x += roots[i] * (long double)kmer[i]; }
        auto a = std::arg(x);
        return pi - 2 * pi / k < a;
    }

    static bool is_decycling_arg_neg(char const* kmer, const uint64_t k) {
        std::complex<long double> x = 0;
        for (uint64_t i = 0; i != k; ++i) { x += roots[i] * (long double)kmer[i]; }
        auto a = std::arg(x);
        return -2 * pi / k < a and a <= 0;
    }
};

template <typename Hasher>
struct double_decycling {
    double_decycling(parameters const& params)
        : m_params(params), m_enum_kmers(m_params.w, m_params.k, m_params.seed) {
        sines.reserve(m_params.k + 1);
        for (uint64_t i = 0; i != m_params.k; ++i) {
            sines.push_back(std::sin(2 * pi * i / m_params.k));
        }
        roots.reserve(m_params.k + 1);
        for (uint64_t i = 0; i != m_params.k; ++i) {
            roots.push_back(std::exp(std::complex<long double>(0, 2 * pi * i / m_params.k)));
        }
    }

    uint64_t sample(char const* window) const {
        using T = std::pair<uint8_t, typename Hasher::hash_type>;
        uint64_t p = -1;
        T min_hash{-1, -1};
        for (uint64_t i = 0; i != m_params.w; ++i) {
            char const* kmer = window + i;
            auto hash = double_decycling_hasher<Hasher>::hash(kmer, m_params.k, m_params.seed);
            if (hash < min_hash) {
                min_hash = hash;
                p = i;
            }
        }
        assert(p < m_params.w);
        return p;
    }

    uint64_t sample(char const* window, bool clear) {
        m_enum_kmers.eat(window, clear);
        return m_enum_kmers.next();
    }

private:
    parameters m_params;
    enumerator<double_decycling_hasher<Hasher>> m_enum_kmers;
};

}  // namespace minimizers