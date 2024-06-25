/* SubModST solver, that solves the Cardinality-Constrained Submodular Monotone
   Subset Maximization problem.
   Copyright (C) 2024  Henning Woydt

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
==============================================================================*/

#ifndef SUBMODST_HEAP_H
#define SUBMODST_HEAP_H

#include <cstdint>
#include <vector>
#include <algorithm>

namespace SubModST {

/**
     * Represents an element on the heap.
     *
     * @tparam SFType Return data type of the score function.
     */
    template<typename SFType>
    struct HeapEntry {
        uint32_t c;
        SFType gain;
    };

    template<typename T>
    std::ostream &operator<<(std::ostream &os,
                             const HeapEntry<T> &entry) {
        os << "(" << entry.c << ", " << entry.gain << ")";
        return os;
    }

    /**
     * Manages the heap for faster pruning.
     *
     * @tparam SFType Return data type of the score function.
     */
    template<typename SFType>
    class Heap {
    private:
        std::vector<HeapEntry<SFType>> m_heap;
        size_t m_heap_max_size;
        size_t m_heap_size;
        SFType m_heap_sum;
        bool abort_early;

    public:
        /**
         * Default constructor.
         */
        Heap() {
            m_heap_max_size = 0;
            m_heap_size = 0;
            m_heap_sum = 0;
            abort_early = false;
        };

        /**
         * Initializes the heap, such that it can hold max_size elements.
         *
         * @param max_size Size of the heap.
         */
        inline void initialize(size_t max_size) {
            m_heap.resize(max_size);
            m_heap_max_size = max_size;
            m_heap_size = 0;
            m_heap_sum = 0;
            abort_early = false;
        }

        /**
         * Pushes the candidate and the score improvement onto the heap.
         *
         * @param c The candidate.
         * @param si The score improvement.
         */
        inline void push(HeapEntry<SFType> e) {
            if (m_heap_size < m_heap_max_size) {
                // enough space available on the heap so push
                m_heap_sum += e.gain;
                m_heap[m_heap_size] = e;
                m_heap_size += 1;

                if (m_heap_size == m_heap_max_size) {
                    // if the heap is full, then heapify
                    std::make_heap(m_heap.begin(), m_heap.end(), [](const HeapEntry<SFType> &a, const HeapEntry<SFType> &b) { return a.gain > b.gain; });
                }
            } else {
                // not enough space available on the heap
                if (e.gain > m_heap[0].gain) {
                    // push to heap is gain is greater
                    m_heap_sum = m_heap_sum - m_heap[0].gain + e.gain;
                    m_heap[0] = e;
                    correct_heap(0);
                }
            }
        }

        /**
         * Updates the score improvement of the entry in the heap.
         *
         * @param e The entry.
         */
        inline void update(HeapEntry<SFType> e) {
            size_t idx = 0;
            for (size_t i = 0; i < m_heap_size; ++i) {
                if (m_heap[i].c == e.c) {
                    m_heap_sum = m_heap_sum - m_heap[i].gain + e.gain;
                    m_heap[i].gain = e.gain;
                    idx = i;
                    break;
                }
            }
            correct_heap(idx);
        }

        /**
         * Edits the heap after changing the gain at index i.
         *
         * @param idx The index that was modified.
         */
        inline void correct_heap(size_t idx) {

            while (true) {
                size_t parent = idx == 0 ? 0 : (idx - 1) / 2;
                size_t left = 2 * idx + 1;
                size_t right = 2 * idx + 2;


                if (m_heap[idx].gain < m_heap[parent].gain) {
                    // swap with parent
                    std::swap(m_heap[idx], m_heap[parent]);
                    idx = parent;
                } else {
                    // no swap with parent, check children
                    size_t smallest = idx;
                    if (left < m_heap_size && m_heap[left].gain < m_heap[smallest].gain) {
                        // swap with left child
                        smallest = left;
                    }
                    if (right < m_heap_size && m_heap[right].gain < m_heap[smallest].gain) {
                        // swap with right child
                        smallest = right;
                    }

                    if (smallest != idx) {
                        std::swap(m_heap[idx], m_heap[smallest]);
                        idx = smallest;
                    } else {
                        break; // heap property is satisfied
                    }
                }
            }
        }

        /**
         * Returns the current sum of the heap.
         *
         * @return The sum.
         */
        inline SFType get_heap_sum() {
            return m_heap_sum;
        };

        /**
         * Returns the minimum score improvement in the heap.
         *
         * @return The minimum score improvement.
         */
        inline SFType get_heap_min() {
            return m_heap[0].gain;
        }

        /**
         * Marks the abort early flag as true.
         */
        inline void mark_abort_early() {
            abort_early = true;
        }

        /**
         * Returns the abort early flag.
         *
         * @return The flag.
         */
        inline bool is_abort_early() {
            return abort_early;
        }
    };
}

#endif //SUBMODST_HEAP_H
