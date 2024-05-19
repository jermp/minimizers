#pragma once

#include <vector>

namespace minimizers {

template <typename T>
struct fixed_size_deque {
    fixed_size_deque() {}
    fixed_size_deque(uint64_t size) : m_begin(0), m_end(0), m_size(0), m_buffer(size) {
        assert(size > 0);
    }

    bool empty() const { return m_size == 0; }

    void push_back(T const& val) {
        assert(m_size < m_buffer.size());
        m_buffer[m_end] = val;
        if (++m_end == m_buffer.size()) m_end = 0;
        ++m_size;
    }

    void pop_front() {
        assert(!empty());
        if (++m_begin == m_buffer.size()) m_begin = 0;
        --m_size;
    }

    void pop_back() {
        assert(!empty());
        if (m_end == 0) m_end = m_buffer.size();
        --m_end;
        --m_size;
    }

    T const& front() const { return m_buffer[m_begin]; }
    T const& back() const { return m_end == 0 ? m_buffer.back() : m_buffer[m_end - 1]; }

private:
    uint64_t m_begin;
    uint64_t m_end;
    uint64_t m_size;
    std::vector<T> m_buffer;
};

template <typename Hasher>
struct enumerator {
    typedef typename Hasher::hash_type hash_type;

    enumerator() {}
    enumerator(uint64_t w, uint64_t k, uint64_t seed)
        : m_w(w), m_k(k), m_seed(seed), m_position(0), m_window(0), m_q(w) {}

    /* return the position of the (leftmost) smallest kmer in the window */
    void eat(char const* window, bool clear) {
        for (uint64_t i = clear ? 0 : m_w - 1; i != m_w; ++i) eat(window + i);
    }

    uint64_t next() {
        uint64_t p = m_q.front().position - m_window;
        assert(p >= 0 and p < m_w);
        m_window += 1;
        return p;
    }

    void eat(char const* kmer) {
        hash_type hash = Hasher::hash(kmer, m_k, m_seed);

        /* Removes from front elements which are no longer in the window */
        while (!m_q.empty() and m_position >= m_w and m_q.front().position <= m_position - m_w) {
            m_q.pop_front();
        }
        /* Removes from back elements which are no longer useful */
        while (!m_q.empty() and hash < m_q.back().hash) m_q.pop_back();

        m_q.push_back({hash, m_position});
        m_position += 1;
    }

    void skip() {
        while (!m_q.empty() and m_position >= m_w and m_q.front().position <= m_position - m_w) {
            m_q.pop_front();
        }
        m_position += 1;
    }

private:
    uint64_t m_w;
    uint64_t m_k;
    uint64_t m_seed;
    uint64_t m_position;
    uint64_t m_window;

    struct kmer_t {
        kmer_t() {}
        kmer_t(hash_type h, uint64_t p) : hash(h), position(p) {}
        hash_type hash;
        uint64_t position;
    };

    fixed_size_deque<kmer_t> m_q;
};

}  // namespace minimizers