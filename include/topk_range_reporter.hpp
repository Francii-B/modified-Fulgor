#pragma once

#include <vector>
#include <algorithm>
#include <cassert>

namespace fulgor {

struct doc_id_range {
    doc_id_range(uint32_t doc_id, uint32_t begin, uint32_t end)
        : doc_id(doc_id), begin(begin), end(end) {}
    doc_id_range() : doc_id_range(-1, -1, -1) {}

    uint32_t doc_id;
    uint32_t begin, end;  // [begin, end]

    bool operator>(doc_id_range const& other) const {
        return end - begin + 1 > other.end - other.begin + 1;
    }
};

struct comparator {
    bool operator()(doc_id_range const& l, doc_id_range const& r) const {
        if (l.end - l.begin == r.end - r.begin) return l.doc_id <= r.doc_id;
        return l > r;
    }
};

struct topk_queue {
    topk_queue(uint32_t k) : m_k(k) {}

    bool insert(doc_id_range const& r) {
        if (m_q.size() < m_k) {
            m_q.push_back(r);
            std::push_heap(m_q.begin(), m_q.end(), comparator());
            return true;
        } else {
            if (r > m_q.front()) {
                std::pop_heap(m_q.begin(), m_q.end(), comparator());
                m_q.back() = r;
                std::push_heap(m_q.begin(), m_q.end(), comparator());
                return true;
            }
        }
        return false;
    }

    void finalize() { std::sort_heap(m_q.begin(), m_q.end(), comparator()); }

    std::vector<doc_id_range> const& topk() const { return m_q; }
    void clear() { m_q.clear(); }

private:
    uint32_t m_k;
    std::vector<doc_id_range> m_q;
};

struct topk_range_reporter {
    static const uint32_t invalid_id = uint32_t(-1);

    struct range {
        range() : begin(invalid_id), end(invalid_id) {}
        uint32_t begin, end;  // [begin, end]
    };

    topk_range_reporter(uint32_t k, uint32_t num_chunks, uint32_t num_docs) {
        m_topk_queues.resize(num_chunks, topk_queue(k));
        m_ranges.resize(num_docs);
    }

    void process_doc_ids(uint32_t chunk_id, std::vector<uint32_t> const& doc_ids) {
        for (auto doc_id : doc_ids) {
            assert(doc_id < m_ranges.size());

            if (m_ranges[doc_id].begin == invalid_id) {  // init range
                assert(m_ranges[doc_id].end == invalid_id);
                m_ranges[doc_id].begin = chunk_id;
                m_ranges[doc_id].end = chunk_id;
                continue;
            }

            if (chunk_id == m_ranges[doc_id].end + 1) {  // extend range
                m_ranges[doc_id].end += 1;
            } else {
                assert(chunk_id > m_ranges[doc_id].end + 1);
                /* process current range and open a new one*/
                process_range(doc_id, m_ranges[doc_id]);
                m_ranges[doc_id].begin = chunk_id;
                m_ranges[doc_id].end = chunk_id;
            }
        }
    }

    void finalize() {
        for (uint64_t doc_id = 0; doc_id != m_ranges.size(); ++doc_id) {
            if (m_ranges[doc_id].begin == invalid_id) continue;
            process_range(doc_id, m_ranges[doc_id]);
        }
        for (auto& q : m_topk_queues) q.finalize();
    }

    std::vector<doc_id_range> const& topk(uint32_t chunk_id) const {
        assert(chunk_id < m_topk_queues.size());
        return m_topk_queues[chunk_id].topk();
    }

    void clear() {
        for (auto& q : m_topk_queues) q.clear();
        for (auto& r : m_ranges) r = range();
    }

private:
    std::vector<topk_queue> m_topk_queues;
    std::vector<range> m_ranges;

    void process_range(uint32_t doc_id, range const& r) {
        for (uint32_t chunk_id = r.begin; chunk_id <= r.end; ++chunk_id) {
            m_topk_queues[chunk_id].insert(doc_id_range(doc_id, r.begin, r.end));
        }
    }
};

}  // namespace fulgor