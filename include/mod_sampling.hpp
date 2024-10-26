#pragma once

#include "external/fastmod/fastmod.h"

namespace minimizers {

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
