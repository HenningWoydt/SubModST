/* CCSMSM solver, that solves the Cardinality-Constrained Submodular Monotone
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

#ifndef CCSMSM_CACHE_H
#define CCSMSM_CACHE_H

#include <map>

namespace CCSMSM {

    /**
     * Represents an entry.
     *
     * @tparam TypeSF Return data type of the score function.
     */
    template<typename SFType>
    struct Pair {
        size_t a;
        size_t b;
        SFType gain;
    };

    template<typename T>
    std::ostream &operator<<(std::ostream &os,
                             const Pair<T> &pair) {
        os << "(" << pair.a << ", " << pair.b << ", " << pair.gain << ")";
        return os;
    }

    /**
     * Returns upper bounds based on 2D marginal gains.
     *
     * @tparam SFType Return data type of the score function.
     */
    template<typename SFType>
    class Cache {
    private:
        std::map<std::pair<size_t, size_t>, SFType> m_cache;

    public:
        /**
         * Default constructor.
         */
        explicit Cache() {
            m_cache.clear();
        }

        /**
         * Clears the cache.
         */
        inline void clear() {
            m_cache.clear();
        }

        /**
         * Inserts a pair into the cache.
         *
         * @param pair The pair.
         */
        inline void insert(const Pair<SFType> &pair) {
            m_cache[std::make_pair(pair.a, pair.b)] = pair.gain;
        }

        /**
         * Retrieves a pair from the cache. Returns true if it exists, and false
         * else. If the pair exists, then pair will hold the gain.
         *
         * @param pair The pair to retrieve.
         * @return True if the pair exists, false else.
         */
        inline bool retrieve(Pair<SFType> &pair) const {
            auto it = m_cache.find(std::make_pair(pair.a, pair.b));
            if (it != m_cache.end()) {
                pair.gain = it->second;
                return true;
            }
            return false;
        }
    };
}

#endif //CCSMSM_CACHE_H
