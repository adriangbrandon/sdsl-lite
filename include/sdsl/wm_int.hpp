/* sdsl - succinct data structures library
    Copyright (C) 2014 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file wm_int.hpp
    \brief wm_int.hpp contains a specialized class for a wavelet tree for
           sequences over large alphabets.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_WM_INT
#define INCLUDED_SDSL_WM_INT

#include "sdsl_concepts.hpp"
#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "select_support_mcl.hpp"
#include "wt_helper.hpp"
#include "util.hpp"
#include <set> // for calculating the alphabet size
#include <map> // for mapping a symbol to its lexicographical index
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <queue>
#include <utility>

#include "init_array.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A wavelet tree class for integer sequences.
/*!
 * \tparam t_bitvector   Type of the bitvector used for representing the wavelet tree.
 * \tparam t_rank        Type of the support structure for rank on pattern `1`.
 * \tparam t_select      Type of the support structure for select on pattern `1`.
 * \tparam t_select_zero Type of the support structure for select on pattern `0`.
 *
 * This wavelet tree variant does not store the two children of a node v aligned
 * with v; it is also known as wavelet matrix.
 *
 * \par References
 *      [1] F. Claude, G. Navarro: ,,The Wavelet Matrix'', Proceedings of
 *          SPIRE 2012.
 *
 *   @ingroup wt
 */
template<class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type>
class wm_int
{
    public:

        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef typename t_bitvector::difference_type difference_type;
        typedef random_access_const_iterator<wm_int> const_iterator;
        typedef const_iterator                       iterator;
        typedef t_bitvector                          bit_vector_type;
        typedef t_rank                               rank_1_type;
        typedef t_select                             select_1_type;
        typedef t_select_zero                        select_0_type;
        typedef wt_tag                               index_category;
        typedef int_alphabet_tag                     alphabet_category;
        enum 	{lex_ordered=0};

        typedef std::pair<value_type, size_type>     point_type;
        typedef std::vector<point_type>              point_vec_type;
        typedef std::pair<size_type, point_vec_type> r2d_res_type;

        struct node_type;


    protected:

        size_type              m_size  = 0;
        size_type              m_sigma = 0;    //<- \f$ |\Sigma| \f$
        bit_vector_type        m_tree;         // bit vector to store the wavelet tree
        rank_1_type            m_tree_rank;    // rank support for the wavelet tree bit vector
        select_1_type          m_tree_select1; // select support for the wavelet tree bit vector
        select_0_type          m_tree_select0;
        uint32_t               m_max_level = 0;
        int_vector<64>         m_zero_cnt;     // m_zero_cnt[i] contains the number of zeros in level i
        int_vector<64>         m_rank_level;   // m_rank_level[i] contains m_tree_rank(i*size())

        void copy(const wm_int& wt)
        {
            m_size          = wt.m_size;
            m_sigma         = wt.m_sigma;
            m_tree          = wt.m_tree;
            m_tree_rank     = wt.m_tree_rank;
            m_tree_rank.set_vector(&m_tree);
            m_tree_select1  = wt.m_tree_select1;
            m_tree_select1.set_vector(&m_tree);
            m_tree_select0  = wt.m_tree_select0;
            m_tree_select0.set_vector(&m_tree);
            m_max_level     = wt.m_max_level;
            m_zero_cnt      = wt.m_zero_cnt;
            m_rank_level    = wt.m_rank_level;
        }

    public:

        const size_type&       sigma = m_sigma;         //!< Effective alphabet size of the wavelet tree.
        const bit_vector_type& tree  = m_tree;          //!< A concatenation of all bit vectors of the wavelet tree.
        const uint32_t&        max_level = m_max_level; //!< Maximal level of the wavelet tree.

        //! Default constructor
        wm_int()
        {
        };

        //! Semi-external constructor
        /*! \param buf         File buffer of the int_vector for which the wm_int should be build.
         *  \param size        Size of the prefix of v, which should be indexed.
         *  \param max_level   Maximal level of the wavelet tree. If set to 0, determined automatically.
         *    \par Time complexity
         *        \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         *        I.e. we need \Order{n\log n} if rac is a permutation of 0..n-1.
         *    \par Space complexity
         *        \f$ n\log|\Sigma| + O(1)\f$ bits, where \f$n=size\f$.
         */
        template<uint8_t int_width>
        wm_int(int_vector_buffer<int_width>& buf, size_type size,
               uint32_t max_level=0) : m_size(size)
        {
            if (0 == m_size)
                return;
            size_type n = buf.size();  // set n
            if (n < m_size) {
                throw std::logic_error("n="+util::to_string(n)+" < "+util::to_string(m_size)+"=m_size");
                return;
            }
            m_sigma = 0; // init sigma

            int_vector<int_width> rac(m_size, 0, buf.width());  // initialize rac

            value_type x = 1;  // variable for the biggest value in rac
            for (size_type i=0; i < m_size; ++i) { // detect the largest value in rac
                if (buf[i] > x)
                    x = buf[i];
                rac[i] = buf[i];
            }

            if (max_level == 0) {
                m_max_level = bits::hi(x)+1; // we need max_level bits to represent all values in the range [0..x]
            } else {
                m_max_level = max_level;
            }


            std::string tree_out_buf_file_name = tmp_file(buf.filename(), "_m_tree");
            osfstream tree_out_buf(tree_out_buf_file_name, std::ios::binary | std::ios::trunc | std::ios::out);   // open buffer for tree
            size_type bit_size = m_size*m_max_level;
            tree_out_buf.write((char*) &bit_size, sizeof(bit_size));    // write size of bit_vector

            std::string zero_buf_file_name = tmp_file(buf.filename(), "_zero_buf");

            size_type tree_pos = 0;
            uint64_t tree_word = 0;

            m_zero_cnt = int_vector<64>(m_max_level, 0); // zeros at level i

            for (uint32_t k=0; k<m_max_level; ++k) {
                uint8_t        width = m_max_level-k-1;
                const uint64_t mask  = 1ULL<<width;
                uint64_t       x     = 0;
                size_type      zeros = 0;
                int_vector_buffer<> zero_buf(zero_buf_file_name, std::ios::out, 1024*1024, m_max_level);
                for (size_t i=0; i<m_size; ++i) {
                    x = rac[i];
                    if (x&mask) {
                        tree_word |= (1ULL << (tree_pos&0x3FULL));
                        zero_buf.push_back(x);
                    } else {
                        rac[zeros++ ] = x;
                    }
                    ++tree_pos;
                    if ((tree_pos & 0x3FULL) == 0) { // if tree_pos % 64 == 0 write old word
                        tree_out_buf.write((char*) &tree_word, sizeof(tree_word));
                        tree_word = 0;
                    }
                }
                m_zero_cnt[k] = zeros;
                for (size_t i=zeros; i<m_size; ++i) {
                    rac[i] = zero_buf[i-zeros];
                }
            }
            if ((tree_pos & 0x3FULL) != 0) { // if tree_pos % 64 > 0 => there are remaining entries we have to write
                tree_out_buf.write((char*) &tree_word, sizeof(tree_word));
            }
            sdsl::remove(zero_buf_file_name);
            tree_out_buf.close();
            m_sigma = std::unique(rac.begin(), rac.end()) - rac.begin();
            rac.resize(0);
            bit_vector tree;
            load_from_file(tree, tree_out_buf_file_name);
            sdsl::remove(tree_out_buf_file_name);
            m_tree = bit_vector_type(std::move(tree));
            util::init_support(m_tree_rank, &m_tree);
            util::init_support(m_tree_select0, &m_tree);
            util::init_support(m_tree_select1, &m_tree);
            m_rank_level = int_vector<64>(m_max_level, 0);
            for (uint32_t k=0; k<m_rank_level.size(); ++k) {
                m_rank_level[k] = m_tree_rank(k*m_size);
            }
        }

        //! Copy constructor
        wm_int(const wm_int& wt)
        {
            copy(wt);
        }

        //! Copy constructor
        wm_int(wm_int&& wt)
        {
            *this = std::move(wt);
        }

        //! Assignment operator
        wm_int& operator=(const wm_int& wt)
        {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Assignment move operator
        wm_int& operator=(wm_int&& wt)
        {
            if (this != &wt) {
                m_size          = wt.m_size;
                m_sigma         = wt.m_sigma;
                m_tree          = std::move(wt.m_tree);
                m_tree_rank     = std::move(wt.m_tree_rank);
                m_tree_rank.set_vector(&m_tree);
                m_tree_select1  = std::move(wt.m_tree_select1);
                m_tree_select1.set_vector(&m_tree);
                m_tree_select0  = std::move(wt.m_tree_select0);
                m_tree_select0.set_vector(&m_tree);
                m_max_level     = std::move(wt.m_max_level);
                m_zero_cnt      = std::move(wt.m_zero_cnt);
                m_rank_level    = std::move(wt.m_rank_level);
            }
            return *this;
        }

        //! Swap operator
        void swap(wm_int& wt)
        {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                std::swap(m_sigma,  wt.m_sigma);
                m_tree.swap(wt.m_tree);
                util::swap_support(m_tree_rank, wt.m_tree_rank, &m_tree, &(wt.m_tree));
                util::swap_support(m_tree_select1, wt.m_tree_select1, &m_tree, &(wt.m_tree));
                util::swap_support(m_tree_select0, wt.m_tree_select0, &m_tree, &(wt.m_tree));
                std::swap(m_max_level,  wt.m_max_level);
                m_zero_cnt.swap(wt.m_zero_cnt);
                m_rank_level.swap(wt.m_rank_level);
            }
        }

        //! Returns the size of the original vector.
        size_type size()const
        {
            return m_size;
        }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const
        {
            return m_size == 0;
        }

        //! Recovers the i-th symbol of the original vector.
        /*! \param i The index of the symbol in the original vector.
         *  \returns The i-th symbol of the original vector.
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        value_type operator[](size_type i)const
        {
            assert(i < size());
            value_type res = 0;
            for (uint32_t k=0; k < m_max_level; ++k) {
                res <<= 1;
                size_type rank_ones = m_tree_rank(i) - m_rank_level[k];
                if (m_tree[i]) { // one at position i => follow right child
                    i = (k+1)*m_size + m_zero_cnt[k] + rank_ones;
                    res |= 1;
                } else { // zero at position i => follow left child
                    auto rank_zeros = (i - k*m_size) - rank_ones;
                    i = (k+1)*m_size + rank_zeros;
                }
            }
            return res;
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *       \f$ \Order{\log |\Sigma|} \f$
         *  \par Precondition
         *       \f$ i \leq size() \f$
         */
        size_type rank(size_type i, value_type c)const
        {
            assert(i <= size());
            if (((1ULL)<<(m_max_level))<=c) { // c is greater than any symbol in wt
                return 0;
            }
            size_type b = 0; // start position of the interval
            uint64_t mask = (1ULL) << (m_max_level-1);
            for (uint32_t k=0; k < m_max_level and i; ++k) {
                size_type rank_b = m_tree_rank(b);
                size_type ones   = m_tree_rank(b + i) - rank_b; // ones in [b..i)
                size_type ones_p = rank_b - m_rank_level[k];    // ones in [level_b..b)
                if (c & mask) { // search for a one at this level
                    i = ones;
                    b = (k+1)*m_size + m_zero_cnt[k] + ones_p;
                } else { // search for a zero at this level
                    i = i-ones;
                    b = (k+1)*m_size + (b - k*m_size - ones_p);
                }
                mask >>= 1;
            }
            return i;
        };

        //! Calculates how many occurrences of symbol wt[i] are in the prefix [0..i-1] of the original sequence.
        /*!
         *  \param i The index of the symbol.
         *  \return  Pair (rank(wt[i],i),wt[i])
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        std::pair<size_type, value_type>
        inverse_select(size_type i)const
        {
            assert(i < size());
            value_type c = 0;
            size_type b = 0; // start position of the interval
            uint64_t mask = (1ULL) << (m_max_level-1);
            for (uint32_t k=0; k < m_max_level; ++k) {
                size_type rank_b = m_tree_rank(b);
                size_type ones   = m_tree_rank(b + i) - rank_b; // ones in [b..i)
                size_type ones_p = rank_b - m_rank_level[k];    // ones in [level_b..b)
                c<<=1;
                if (m_tree[b+i]) { // go to the right child
                    i = ones;
                    b = (k+1)*m_size + m_zero_cnt[k] + ones_p;
                    c|=1;
                } else { // go to the left child
                    i = i-ones;
                    b = (k+1)*m_size + (b - k*m_size - ones_p);
                }
                mask >>= 1;
            }
            return std::make_pair(i,c);
        }

        //! Calculates the i-th occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The i-th occurrence.
         *  \param c The symbol c.
         *  \par Time complexity
         *       \f$ \Order{\log |\Sigma|} \f$
         *  \par Precondition
         *       \f$ 1 \leq i \leq rank(size(), c) \f$
         */
        size_type select(size_type i, value_type c)const
        {
            //assert(1 <= i and i <= rank(size(), c));
            uint64_t mask = 1ULL << (m_max_level-1);
            int_vector<64> m_path_off(max_level+1);
            int_vector<64> m_path_rank_off(max_level+1);
            m_path_off[0] = m_path_rank_off[0] = 0;
            size_type b = 0; // start position of the interval
            size_type r = i;
            for (uint32_t k=0; k < m_max_level and i; ++k) {
                size_type rank_b = m_tree_rank(b);
                size_type ones   = m_tree_rank(b + r) - rank_b; // ones in [b..i)
                size_type ones_p = rank_b - m_rank_level[k];    // ones in [0..b)
                if (c & mask) { // search for a one at this level
                    r = ones;
                    b = (k+1)*m_size + m_zero_cnt[k] + ones_p;
                } else { // search for a zero at this level
                    r = r-ones;
                    b = (k+1)*m_size + (b - k*m_size - ones_p);
                }
                mask >>= 1;
                m_path_off[k+1] = b;
                m_path_rank_off[k] = rank_b;
            }
            mask = 1ULL;
            for (uint32_t k=m_max_level; k>0; --k) {
                b = m_path_off[k-1];
                size_type rank_b = m_path_rank_off[k-1];
                if (c & mask) { // right child => search i'th one
                    i = m_tree_select1(rank_b + i) - b + 1;
                } else { // left child => search i'th zero
                    i = m_tree_select0(b - rank_b + i) - b + 1;
                }
                mask <<= 1;
            }
            return i-1;
        };

	std::pair<size_type, size_type> select_next(size_type i, value_type c, size_type n_elems) const
        {
            //assert(1 <= i and i <= rank(size(), c));
            assert(i <= size());
	        uint64_t mask = 1ULL << (m_max_level-1);
            int_vector<64> m_path_off(max_level+1);
            int_vector<64> m_path_rank_off(max_level+1);
            m_path_off[0] = m_path_rank_off[0] = 0;
            size_type b = 0; // start position of the interval
            size_type r = i;
            for (uint32_t k=0; k < m_max_level and i; ++k) {
                size_type rank_b = m_tree_rank(b);
                size_type ones   = m_tree_rank(b + r) - rank_b; // ones in [b..i)
                size_type ones_p = rank_b - m_rank_level[k];    // ones in [0..b)
                if (c & mask) { // search for a one at this level
                    r = ones;
                    b = (k+1)*m_size + m_zero_cnt[k] + ones_p;
                } else { // search for a zero at this level
                    r = r-ones;
                    b = (k+1)*m_size + (b - k*m_size - ones_p);
                }
                mask >>= 1;
                m_path_off[k+1] = b;
                m_path_rank_off[k] = rank_b;
            }
            mask = 1ULL;
	    i = r+1;
	    if (i > n_elems) return std::pair<size_type, size_type>(0,0);
            for (uint32_t k=m_max_level; k>0; --k) {
                b = m_path_off[k-1];
                size_type rank_b = m_path_rank_off[k-1];
                if (c & mask) { // right child => search i'th one
                    i = m_tree_select1(rank_b + i) - b + 1;
                } else { // left child => search i'th zero
                    i = m_tree_select0(b - rank_b + i) - b + 1;
                }
                mask <<= 1;
            }
            return std::pair<size_type, size_type>(i-1, r);
        };

    uint64_t select_next_ranges_node(node_type &v, const std::vector<range_type> &ranges, range_type sigma_range, value_type c) {
        /*if (!overlaps(sigma_range, sigma_ranges)) {
            return -1ULL;
        }*/

        if (is_leaf(v)) return v.offset + ranges[0][0];

        std::vector<range_type> left_ranges, right_ranges;
        auto child = my_expand_ranges(v, ranges, left_ranges, right_ranges);
        size_type mid = (sigma_range[0] + sigma_range[1]+1)>>1;
        uint64_t p = -1ULL;
        if (c < mid && !left_ranges.empty()) {
            p = select_next_ranges_node(child[0], left_ranges, {sigma_range[0], mid-1}, c);
            if (p != -1ULL) {
                p = p - (v.level+1) * m_size + 1;  //corresponds with the p-th 0-bit of v.level
                p = m_tree_select0((v.level*m_size - m_rank_level[v.level])+p); //position in the bitvector (v.level)
            }
        }else if (c >= mid && !right_ranges.empty()) {
            p = select_next_ranges_node(child[1], right_ranges, {mid, sigma_range[1]}, c);
            if (p != -1ULL) {
                p = p - (v.level+1) * m_size - m_zero_cnt[v.level] +1; //corresponds with the p-th 1-bit of v.level
                p = m_tree_select1(m_rank_level[v.level]+p); //position in the bitvector (v.level)
            }
        }
        return p;
    }

    /***
    * The position of ranges whose contains the first symbol c
    */
    size_type select_next_ranges(const std::vector<range_type> &ranges, value_type c) {
        node_type v = root();
        range_type sigma_range = {0, (1ULL<<m_max_level)-1};
        auto i = select_next_ranges_node(v, ranges, sigma_range, c);
        if (i == -1ULL) return 0;
        return i;
    }



    //Implemented by Adrian Gomez-Brandon

    bool overlaps(const range_type &range, const std::vector<range_type> &sigma_ranges) {

        for (const auto &r : sigma_ranges) {
            if (r[0] <= range[1] && r[1] >= range[0]) {
                return true;
            }
        }
        return false;
    }

    //Computes the state of the range in sigma_ranges
    //- 0: not overlapping
    //- 1: overlapping
    //- 2: contained
    uint check_range(const range_type &range, const std::vector<range_type> &sigma_ranges) {
        for (const auto &r : sigma_ranges) {
            if (r[0] <= range[0] && r[1] >= range[1]) {
                return 2;
            }
            if (r[0] <= range[1] && r[1] >= range[0]) {
                return 1;
            }
        }
        return 0;
    }

    uint64_t select_next_node(node_type &v, range_type range, range_type sigma_range,  const std::vector<range_type> &sigma_ranges) {

        /*if (!overlaps(sigma_range, sigma_ranges)) {
            return -1ULL;
        }*/

        if (is_leaf(v)) return v.offset + range[0];

        uint64_t rank_b;
        range_type left_range, right_range;
        auto child = my_expand(v, range, left_range, right_range, rank_b);
        uint64_t lp = -1ULL, rp = -1ULL;
        size_type mid = (sigma_range[0] + sigma_range[1]+1)>>1;
        uint lstate, rstate;

        if (!sdsl::empty(left_range)) {
            lstate = check_range({sigma_range[0], mid-1}, sigma_ranges);
            if (lstate == 2) {
                lp = child[0].offset + left_range[0]; //position in the bitvector (v.level+1)
                lp = lp - (v.level+1) * m_size + 1;  //corresponds with the lp-th 0-bit of v.level
                lp = m_tree_select0((v.level*m_size - m_rank_level[v.level])+lp); //position in the bitvector (v.level)
            }else if (lstate == 1) {
                lp = select_next_node(child[0], left_range, {sigma_range[0], mid-1}, sigma_ranges); //position in the bitvector (v.level+1)
                if (lp != -1ULL) {
                    lp = lp - (v.level+1) * m_size + 1; //corresponds with the lp-th 0-bit of v.level
                    lp = m_tree_select0((v.level*m_size - m_rank_level[v.level])+lp); //position in the bitvector (v.level)
                }
            }
        }
        if (!sdsl::empty(right_range)) {
            rstate = check_range({mid, sigma_range[1]}, sigma_ranges);
            if (rstate == 2) {
                rp = child[1].offset + right_range[0]; //position in the bitvector (v.level+1)
                rp = rp - (v.level+1) * m_size - m_zero_cnt[v.level] +1; //corresponds with the rp-th 1-bit of v.level
                rp = m_tree_select1(m_rank_level[v.level]+rp); //position in the bitvector (v.level)
            }else if (rstate == 1) {
                rp = select_next_node(child[1], right_range, {mid, sigma_range[1]}, sigma_ranges); //position in the bitvector (v.level+1)
                if (rp != -1ULL) {
                    rp = rp - (v.level+1) * m_size - m_zero_cnt[v.level] +1; //corresponds with the rp-th 1-bit of v.level
                    rp = m_tree_select1(m_rank_level[v.level]+rp); //position in the bitvector (v.level)
                }
            }
        }

        /*if (lp != -1ULL && rp != -1ULL) {
            lp = lp - (v.level+1) * m_size + 1;
            rp = rp - (v.level+1) * m_size - m_zero_cnt[v.level] +1;
            res = std::min(m_tree_select0((v.level*m_size - m_rank_level[v.level])+lp), m_tree_select1(m_rank_level[v.level]+rp));
        } else if (lp != -1ULL) {
           lp = lp - (v.level+1) * m_size + 1;
            res = m_tree_select0((v.level*m_size - m_rank_level[v.level])+lp);
        } else if (rp != -1ULL) {
            rp = rp - (v.level+1) * m_size - m_zero_cnt[v.level] +1;
            res =m_tree_select1(m_rank_level[v.level]+rp);
        }*/
        return std::min(lp, rp);

    }
    /***
     * The position of the range [p, \infty) whose symbol is in any range of sigma_ranges
     */
     size_type  select_next(uint64_t p, std::vector<sdsl::range_type> &sigma_ranges) {
        node_type v = root();
        range_type range = {p, size()-1};
        range_type sigma_range = {0, (1ULL<<m_max_level)-1};
        if (!overlaps(sigma_range, sigma_ranges)) return 0;
        auto i = select_next_node(v, range, sigma_range, sigma_ranges);
        if (i == -1ULL) return 0;
        return i;
    }


    //Implemented by Adrian Gomez-Brandon
    //Given a node v and its i-th element, return the value of the symbol at position i
    value_type value_node(node_type &v, size_type i) {
        value_type res = v.sym;
        i = v.offset + i; //position in the bitvector
        for (uint32_t k=v.level; k < m_max_level; ++k) {
            res <<= 1;
            size_type rank_ones = m_tree_rank(i) - m_rank_level[k];
            if (m_tree[i]) { // one at position i => follow right child
                i = (k+1)*m_size + m_zero_cnt[k] + rank_ones;
                res |= 1;
            } else { // zero at position i => follow left child
                auto rank_zeros = (i - k*m_size) - rank_ones;
                i = (k+1)*m_size + rank_zeros;
            }
        }
        return res;
    }


    //Implemented by Adrian Gomez-Brandon
    std::pair<size_type, value_type>  select_next_pos_with_value_node(node_type &v, range_type range, range_type sigma_range,  const std::vector<range_type> &sigma_ranges) {

        if (is_leaf(v)) return {v.offset + range[0], v.sym};

        uint64_t rank_b;
        range_type left_range, right_range;
        auto child = my_expand(v, range, left_range, right_range, rank_b);
        std::pair<size_type, value_type> lp = {-1ULL, 0}, rp =  {-1ULL, 0};
        size_type mid = (sigma_range[0] + sigma_range[1]+1)>>1;
        uint lstate, rstate;

        if (!sdsl::empty(left_range)) {
            lstate = check_range({sigma_range[0], mid-1}, sigma_ranges);
            if (lstate == 2) { //contained
                lp.first = child[0].offset + left_range[0]; //position in the bitvector (v.level+1)
                //TODO: find value
                lp.second = value_node(child[0], left_range[0]);
                lp.first = lp.first - (v.level+1) * m_size + 1;  //corresponds with the lp-th 0-bit of v.level
                lp.first = m_tree_select0((v.level*m_size - m_rank_level[v.level])+lp.first); //position in the bitvector (v.level)

            }else if (lstate == 1) { //overlaps
                lp = select_next_pos_with_value_node(child[0], left_range, {sigma_range[0], mid-1}, sigma_ranges); //position in the bitvector (v.level+1)
                if (lp.first != -1ULL) {
                    lp.first = lp.first - (v.level+1) * m_size + 1; //corresponds with the lp-th 0-bit of v.level
                    lp.first = m_tree_select0((v.level*m_size - m_rank_level[v.level])+lp.first); //position in the bitvector (v.level)
                }
            }
        }
        if (!sdsl::empty(right_range)) {
            rstate = check_range({mid, sigma_range[1]}, sigma_ranges);
            if (rstate == 2) { //contained
                rp.first = child[1].offset + right_range[0]; //position in the bitvector (v.level+1)
                //TODO: find value
                rp.second = value_node(child[1], right_range[0]);
                rp.first = rp.first - (v.level+1) * m_size - m_zero_cnt[v.level] +1; //corresponds with the rp-th 1-bit of v.level
                rp.first = m_tree_select1(m_rank_level[v.level]+rp.first); //position in the bitvector (v.level)
            }else if (rstate == 1) { //overlaps
                rp = select_next_pos_with_value_node(child[1], right_range, {mid, sigma_range[1]}, sigma_ranges); //position in the bitvector (v.level+1)
                if (rp.first != -1ULL) {
                    rp.first = rp.first - (v.level+1) * m_size - m_zero_cnt[v.level] +1; //corresponds with the rp-th 1-bit of v.level
                    rp.first = m_tree_select1(m_rank_level[v.level]+rp.first); //position in the bitvector (v.level)
                }
            }
        }
        if (lp.first > rp.first) return rp;
        return lp;

    }

    /***
     * The next position and value of the range [p, \infty) whose symbol is in any range of sigma_ranges
     */
    std::pair<size_type, value_type>  select_next_pos_with_value(uint64_t p, std::vector<sdsl::range_type> &sigma_ranges) {
        node_type v = root();
        range_type root_range = {p, size()-1};
        range_type root_sigma_range = {0, (1ULL<<m_max_level)-1};
        if (!overlaps(root_sigma_range, sigma_ranges)) return {0, 0};
        auto pv = next_pos_value_node(v, root_range, root_sigma_range, sigma_ranges);
        if (pv.first == -1ULL) return {0, 0};
        return pv;
    }

    //Implemented by Adrian Gomez-Brandon
    std::pair<value_type, size_type> range_min_value_node(node_type &v, const std::vector<range_type> &ranges, const range_type &sigma_range) {

        if (is_leaf(v)) return {v.sym, v.offset};

        std::vector<range_type> left_ranges, right_ranges;
        auto child = my_expand_ranges(v, ranges, left_ranges, right_ranges);
        size_type mid = (sigma_range[0] + sigma_range[1]+1)>>1;
        if (!left_ranges.empty()) {
            return range_min_value_node(child[0], left_ranges, {sigma_range[0], mid-1});
        }else {
            return range_min_value_node(child[1], right_ranges, {mid, sigma_range[1]});
        }
    }

    std::pair<value_type, size_type> range_next_value_node(node_type &v, const value_type val, const std::vector<range_type> &ranges, const range_type &sigma_range) {

        if (is_leaf(v)) return {v.sym, v.offset};

        std::vector<range_type> left_ranges, right_ranges;
        auto child = my_expand_ranges(v, ranges, left_ranges, right_ranges);
        size_type mid = (sigma_range[0] + sigma_range[1]+1)>>1;
        if (!left_ranges.empty() && val < mid) {
            auto aux =  range_next_value_node(child[0], val, left_ranges, {sigma_range[0], mid-1});
            if (aux) return aux;
            if (!right_ranges.empty()) return range_min_value_node(child[1], right_ranges, {mid, sigma_range[1]});
        }else if (!right_ranges.empty()){
            return range_next_value_node(child[1], val, right_ranges, {mid, sigma_range[1]});
        }
        return {0,0};
    }

    /***
     * The smallest value in [val, \infty) placed in any range of ranges and rank(i) of that value
     */
    std::pair<value_type, size_type> range_next_value(const value_type &val, const std::vector<sdsl::range_type> &ranges) {
        node_type v = root();
        range_type sigma_range = {0, (1ULL<<m_max_level)-1};
        if (val > sigma_range[1]) return {0, 0};
        return range_next_value_node(v, val, ranges, sigma_range);
    }

    /***
    * The smallest value placed in any range of ranges and rank(i) of that value
    */
    std::pair<value_type, size_type> range_min_value(const std::vector<sdsl::range_type> &ranges) {
        node_type v = root();
        range_type sigma_range = {0, (1ULL<<m_max_level)-1};
        return range_min_value_node(v, ranges, sigma_range);
    }


    bool in_sigma_range(const range_type &range, const std::vector<range_type> &sigma_ranges, size_type i_sr) {
        while (sigma_ranges[i_sr][1] < range[0] && i_sr < sigma_ranges.size()) {
            ++i_sr;
        }
        return (i_sr < sigma_ranges.size() && sigma_ranges[i_sr][0] <= range[1] && sigma_ranges[i_sr][1] >= range[0]);
    }

    void range2d_values_node(node_type &v, const sdsl::range_type &range, const std::vector<range_type> &sigma_ranges,
                             const range_type &sigma_range, size_type &i_sr, std::vector<value_type> &res) {

        if (!in_sigma_range(sigma_range, sigma_ranges, i_sr)) {
            return;
        }

        if (is_leaf(v)) {
            res.push_back(v.sym);
        }else{

            size_type rnk;
            range_type left_range, right_range;
            auto child =  my_expand(v, range,
                                                  left_range, right_range, rnk);
            size_type mid = (sigma_range[0] + sigma_range[1]+1)>>1;
            if(!sdsl::empty(left_range)){
                range2d_values_node(child[0], left_range, sigma_ranges, {sigma_range[0], mid-1}, i_sr, res);
            }
            if(!sdsl::empty(right_range)){
                range2d_values_node(child[1], right_range, sigma_ranges, {mid, sigma_range[1]}, i_sr, res);
            }
        }

    }


    //Implemented by Adrian Gomez-Brandon
    //Pre: always the same range of positions and different ranges of symbols
    std::vector<value_type> range2d_values(const sdsl::range_type &range,
                                           const std::vector<sdsl::range_type> &sigma_ranges) {

        std::vector<value_type> result;
        node_type v = root();
        range_type sigma_range = {0, (1ULL<<m_max_level)-1};
        size_type i_sr = 0; //current sigma range to check
        range2d_values_node(v, range, sigma_ranges, sigma_range, i_sr, result);
        return result;

    }



	value_type range_minimum_query(size_type i, size_type j) const
	{
	    return _range_minimum_query(i, j, 0, 0, 0);	
	}

	value_type _range_minimum_query(size_type i, size_type j,
                                         uint32_t depth, size_type b,
                                         value_type res) const
	{
	    if (depth == m_max_level)
	        return res;
	    else {
	        size_type rank_0_b = m_tree_rank(b); // ones in [0..b)
	        size_type rank_b_i = m_tree_rank(b + i) - rank_0_b; // ones in [b..i)
	        size_type rank_b_j = m_tree_rank(b + j + 1) - rank_0_b;
		    size_type ones_p   = rank_0_b - m_rank_level[depth];

            size_type i_l = i - rank_b_i; // zeroes in [b..i)
            size_type j_l = j - rank_b_j;
            size_type i_r = i - i_l;
            size_type j_r = j - 1 - j_l;
            size_type n_l = j_l - i_l + 1;

            res <<= 1;
            if (n_l == 0) { // no left child, recurse on the rigth child
                b = (depth+1)*m_size + m_zero_cnt[depth] + ones_p;
                res |= 1;
                return _range_minimum_query(i_r, j_r, depth+1, b, res);
            }
            else { // recurse on the left child
                b = (depth+1)*m_size + (b - depth*m_size - ones_p);
                return _range_minimum_query(i_l, j_l, depth+1, b, res);
            }
	    }
	}

        // Implemented by Diego Arroyuelo
	value_type range_next_value(value_type x, size_type i, size_type j) const
	{
            if (((1ULL)<<(m_max_level))<=x) { // c is greater than any symbol in wt
                return 0;
            }
            return _range_next_value(x, i, j, 0, 0, 0);
        };


	value_type _range_next_value_min(size_type i, size_type j, 
                                         uint32_t depth, size_type b,
                                         value_type res) const
        {
            if (b+i > b+j) return 0;

            if (depth == m_max_level)
                return res;
            else {
                size_type rank_0_b = m_tree_rank(b); // ones in [0..b)
                size_type rank_b_i = m_tree_rank(b + i) - rank_0_b; // ones in [b..i) 
                size_type rank_b_j = m_tree_rank(b + j + 1) - rank_0_b;
                size_type ones_p   = rank_0_b - m_rank_level[depth];
 
                size_type i_l = i - rank_b_i; // zeroes in [b..i)
                size_type j_l = j - rank_b_j;
                size_type i_r = i - i_l;
                size_type j_r = j - 1 - j_l;
                size_type n_l = j_l - i_l + 1;

                res <<= 1;
                if (n_l == 0) { // no left child, recurse on the rigth child
                    b = (depth+1)*m_size + m_zero_cnt[depth] + ones_p;
                    res |= 1;
                    return _range_next_value_min(i_r, j_r, depth+1, b, res);
                } else { // recurse on the left child
                    b = (depth+1)*m_size + (b - depth*m_size - ones_p);	
                    return _range_next_value_min(i_l, j_l, depth+1, b, res);
                }
            }
        };


        value_type _range_next_value(value_type x, size_type i, size_type j, 
                                     uint32_t depth, size_type b,
                                     value_type res) const
        {
            if (b+i > b+j) 
                return 0;
            else 
                if (depth == m_max_level)
                    return res;
                else {
                    size_type rank_0_b = m_tree_rank(b); // ones in [0..b)
                    size_type rank_b_i = m_tree_rank(b + i) - rank_0_b; // ones in [b..i) 
                    size_type rank_b_j = m_tree_rank(b + j + 1) - rank_0_b;
                    size_type ones_p   = rank_0_b - m_rank_level[depth];
		            size_type i_l = i - rank_b_i; // zeroes in [b..i)
                    size_type j_l = j - rank_b_j;
                    size_type i_r = i - i_l;
                    size_type j_r = j - 1 - j_l;

                    uint64_t mask = (1ULL) << (m_max_level-1-depth);
                    res <<= 1;
                    if (x & mask) { // recurse on the rigth child
                        b = (depth+1)*m_size + m_zero_cnt[depth] + ones_p;
                        res |= 1;
                        return _range_next_value(x, i_r, j_r, depth+1, b, res);
                    } else { // recurse on the left child
                        b = (depth+1)*m_size + (b - depth*m_size - ones_p); 
                        value_type y = _range_next_value(x, i_l, j_l, depth+1, b, res);

                        if (y != 0)
                            return y;
                        else {
                            b = (depth+1)*m_size + m_zero_cnt[depth] + ones_p;
                            res |= 1;
                            return _range_next_value_min(i_r, j_r, depth+1, b, res);
                        }
                    }
                }					
        };


        // Implemented by Diego Arroyuelo
        // Given a value x and a range [i, j], returns  a pair (r, k) such that k is the lefmost index such that i \le k \le j 
        // such that s[k] = r \ge x holds (here s is the indexed sequence). 
        // The second component is j+1 if x is greater than any value in s[i..j].
        std::pair<value_type,size_type> range_next_value_pos(value_type x, size_type i, size_type j) 
        {
            if (((1ULL)<<(m_max_level))<=x) { // c is greater than any symbol in wt
                return std::pair<value_type,size_type>((value_type)-1,(size_type)-1);
            }
            size_type pos;
            value_type temp = _range_next_value_pos(x, i, j, 0, 0, 0, pos);
            return std::pair<value_type,size_type>(temp, pos-1);
        };


        value_type _range_next_value_min_pos(size_type i, size_type j,
                                              uint32_t depth, size_type b,
                                              value_type res, size_type& pos)
        {
            if (b+i > b+j) {
                pos = j+2;
                return 0;
            }

            if (depth == m_max_level) {
	        pos = i+1;
                return res;
            }		
            else {
                size_type rank_0_b = m_tree_rank(b); // ones in [0..b)
                size_type rank_b_i = m_tree_rank(b + i) - rank_0_b; // ones in [b..i)
                size_type rank_b_j = m_tree_rank(b + j + 1) - rank_0_b;
                size_type ones_p   = rank_0_b - m_rank_level[depth];

                size_type i_l = i - rank_b_i; // zeroes in [b..i)
                size_type j_l = j - rank_b_j;
                size_type i_r = i - i_l;
                size_type j_r = j - 1 - j_l;
                size_type n_l = j_l - i_l + 1;

                res <<= 1;
                if (n_l == 0) { // no left child, recurse on the rigth child
                    size_type ant_b = b;
                    b = (depth+1)*m_size + m_zero_cnt[depth] + ones_p;
                    res |= 1;
                    value_type temp = _range_next_value_min_pos(i_r, j_r, depth+1, b, res, pos);
		    pos = m_tree_select1(rank_0_b + pos) - ant_b + 1;
		    return temp;
                } else { // recurse on the left child
                    size_type ant_b = b;
                    b = (depth+1)*m_size + (b - depth*m_size - ones_p);
                    value_type temp = _range_next_value_min_pos(i_l, j_l, depth+1, b, res, pos);
                    pos = m_tree_select0(pos + ant_b-(rank_0_b)) - ant_b + 1;
		    return temp;
                }
            }
        };

        value_type _range_next_value_pos(value_type x, size_type i, size_type j,
                                          uint32_t depth, size_type b,
                                          value_type res, size_type& pos) 
        {
            if (b+i > b+j) {
                pos = j+2;
                return 0;
            }
            else
                if (depth == m_max_level) {
		    pos = i+1;
                    return res;
                } else {
                    size_type rank_0_b = m_tree_rank(b); // ones in [0..b)
                    size_type rank_b_i = m_tree_rank(b + i) - rank_0_b; // ones in [b..i)
                    size_type rank_b_j = m_tree_rank(b + j + 1) - rank_0_b;
                    size_type ones_p   = rank_0_b - m_rank_level[depth];
                    size_type i_l = i - rank_b_i; // zeroes in [b..i)
                    size_type j_l = j - rank_b_j;
                    size_type i_r = i - i_l;
                    size_type j_r = j - 1 - j_l;

                    uint64_t mask = (1ULL) << (m_max_level-1-depth);
                    res <<= 1;
                    if (x & mask) { // recurse on the rigth child
                        size_type ant_b = b;
                        b = (depth+1)*m_size + m_zero_cnt[depth] + ones_p;
                        res |= 1;
                        value_type temp = _range_next_value_pos(x, i_r, j_r, depth+1, b, res, pos);
                        if (temp)
			    pos = m_tree_select1(rank_0_b + pos) - ant_b + 1;
                        else
                            pos = j+2;
                        return temp;
                        
                    } else { // recurse on the left child
                        size_type ant_b = b;
                        b = (depth+1)*m_size + (b - depth*m_size - ones_p);
                        value_type y = _range_next_value_pos(x, i_l, j_l, depth+1, b, res, pos);

                        if (y != 0) {
                            pos = m_tree_select0(pos + ant_b-(rank_0_b)) - ant_b + 1;
                            return y;
                        }
                        else {
                            b = (depth+1)*m_size + m_zero_cnt[depth] + ones_p;
                            res |= 1;
                            value_type temp = _range_next_value_min_pos(i_r, j_r, depth+1, b, res, pos);
                            pos = m_tree_select1(pos + rank_0_b) - ant_b + 1;
			    return temp;
                        }
                    }
                }
        };



        // Implemented by Diego Arroyuelo
        template<typename word_t>
        void mark(value_type c, std::vector<word_t>& B_wt, word_t B_c)const
        {
            //assert(1 <= i and i <= rank(size(), c));
            size_type i = 1; // OJO AQUI
            uint64_t mask = 1ULL << (m_max_level-1);
            int_vector<64> m_path_off(max_level+1);
            int_vector<64> m_path_rank_off(max_level+1);
            m_path_off[0] = m_path_rank_off[0] = 0;
            size_type b = 0; // start position of the interval
            size_type r = i;
            uint64_t pos = 1; // assume positions in B_wt array start at 1
            for (uint32_t k=0; k < m_max_level and i; ++k) {
                size_type rank_b = m_tree_rank(b); 
                size_type ones   = m_tree_rank(b + r) - rank_b; // ones in [b..i)
                size_type ones_p = rank_b - m_rank_level[k];    // ones in [0..b)
                B_wt[pos] |= B_c;
                if (c & mask) { // search for a one at this level
                    r = ones;
                    b = (k+1)*m_size + m_zero_cnt[k] + ones_p;
                    pos = 2*pos+1;
                } else { // search for a zero at this level
                    r = r-ones;
                    b = (k+1)*m_size + (b - k*m_size - ones_p);
                    pos = 2*pos;
                }
                mask >>= 1;
                m_path_off[k+1] = b; 
                m_path_rank_off[k] = rank_b;
            }
            B_wt[pos] |= B_c;   // esto no deberia ser necesario si al buscar uso el arreglo B
            //std::cout << "Mark con pos=" << pos << std::endl;
        };


        // Implemented by Diego Arroyuelo
        template<typename word_t>
        void unmark(value_type c, std::vector<word_t>& B_wt)const
        {
            //assert(1 <= i and i <= rank(size(), c));
            size_type i = 1; // OJO AQUI
            uint64_t mask = 1ULL << (m_max_level-1);
            int_vector<64> m_path_off(max_level+1);
            int_vector<64> m_path_rank_off(max_level+1);
            m_path_off[0] = m_path_rank_off[0] = 0;
            size_type b = 0; // start position of the interval
            size_type r = i;
            uint64_t pos = 1; // assume positions in B_wt array start at 1
            for (uint32_t k=0; k < m_max_level and i; ++k) {
                size_type rank_b = m_tree_rank(b); 
                size_type ones   = m_tree_rank(b + r) - rank_b; // ones in [b..i)
                size_type ones_p = rank_b - m_rank_level[k];    // ones in [0..b)
                B_wt[pos] &= 0;
                if (c & mask) { // search for a one at this level
                    r = ones;
                    b = (k+1)*m_size + m_zero_cnt[k] + ones_p;
                    pos = 2*pos+1;
                } else { // search for a zero at this level
                    r = r-ones;
                    b = (k+1)*m_size + (b - k*m_size - ones_p);
                    pos = 2*pos;
                }
                mask >>= 1;
                m_path_off[k+1] = b;
                m_path_rank_off[k] = rank_b;
            }
        };
        

        // Implemented by Diego Arroyuelo
        std::vector<value_type>
        all_values_in_range(size_type lb, size_type rb, bool report=true) const
        {
            size_type cnt_answers = 0;
            std::vector<value_type> res_vec;
            if (lb <= rb) {
                std::vector<size_type> is(m_max_level+1);
                std::vector<size_type> rank_off(m_max_level+1);
                _all_values_in_range(root(), {lb, rb}, 0, is,
                                 rank_off, res_vec, report, cnt_answers);
            }
            return res_vec;
        }

        // Implemented by Diego Arroyuelo
        void 
        all_values_in_range_bounded(size_type lb, size_type rb, 
                                    std::vector<value_type>& res_vec, 
                                    size_type bound, bool report=true) const
        {
	    size_type cnt_answers = 0;
	    if (lb <= rb)  {
		_all_values_in_range_bounded(root(), {lb, rb}, 0, res_vec, report, 
                                             cnt_answers, bound);
	    }
	}
	
	


        // implemented by Diego Arroyuelo
        void
        _all_values_in_range(node_type v, range_type r,
                         size_type ilb, std::vector<size_type>& is,
                         std::vector<size_type>& rank_off, 
                         std::vector<value_type>& res_vec,
                         bool report, size_type& cnt_answers)
        const
        {
            
            using std::get;
            if (get<0>(r) > get<1>(r))
                return;
            is[v.level] = v.offset + get<0>(r);

            if (v.level == m_max_level) {
                res_vec.emplace_back(v.sym);
                cnt_answers += 1; 
                return;
            } else {
                rank_off[v.level] = m_tree_rank(is[v.level]);
            }
            size_type irb = ilb + (1ULL << (m_max_level-v.level));
            size_type mid = (irb + ilb)>>1;

            auto c_v = expand(v);
            auto c_r = expand(v, r);

            if (!sdsl::empty(get<0>(c_r)) and mid) {
                _all_values_in_range(get<0>(c_v),get<0>(c_r), ilb, is, 
                                 rank_off, res_vec, report, cnt_answers);
            }
            if (!sdsl::empty(get<1>(c_r))) {
                _all_values_in_range(get<1>(c_v), get<1>(c_r), mid, is, 
                                     rank_off, res_vec, report, cnt_answers);
            }
        }


        // implemented by Diego Arroyuelo
        void
        _all_values_in_range_bounded(node_type v, range_type r,
                         size_type ilb, //std::vector<size_type>& is,
                         //std::vector<size_type>& rank_off,
                         std::vector<value_type>& res_vec,
                         bool report, size_type& cnt_answers, size_type bound)
        const
        {

            using std::get;
            if ((get<0>(r) > get<1>(r)) or (cnt_answers >= bound))
                return;

            if (v.level == m_max_level) {
                res_vec.emplace_back(v.sym);
                cnt_answers += 1;
                return;
            } 
              
            size_type irb = ilb + (1ULL << (m_max_level-v.level));
            size_type mid = (irb + ilb)>>1;

            auto c_v = expand(v);
            auto c_r = expand(v, r);

            if (!sdsl::empty(get<0>(c_r)) and mid) {
                _all_values_in_range_bounded(get<0>(c_v),get<0>(c_r), ilb, //is,
                                 /*rank_off,*/ res_vec, report, cnt_answers, bound);
            }
            if (!sdsl::empty(get<1>(c_r))) {
                _all_values_in_range_bounded(get<1>(c_v), get<1>(c_r), mid, /*is,*/
                                     /*rank_off,*/ res_vec, report, cnt_answers, bound);
            }
        }

        // Implemented by Diego Arroyuelo
        template<typename word_t>
        std::vector<value_type>
        all_active_p_values_in_range(size_type lb, size_type rb, 
                                   std::vector<word_t>& B_wt, word_t D, 
                                   bool report=true) const
        {

            size_type cnt_answers = 0;
            std::vector<value_type> res_vec;
            if (lb <= rb) {
                std::vector<size_type> is(m_max_level+1);
                std::vector<size_type> rank_off(m_max_level+1); 
                _all_active_p_values_in_range<word_t>(root(), {lb, rb}, 0, is,
                                 rank_off, res_vec, report, cnt_answers, B_wt, D, 1);
            }
            return res_vec;
        }

        // implemented by Diego Arroyuelo
        template<typename word_t>
        void
        _all_active_p_values_in_range(node_type v, range_type r, size_type ilb, std::vector<size_type>& is,
                         std::vector<size_type>& rank_off, std::vector<value_type>& res_vec,
                         bool report, size_type& cnt_answers, std::vector<word_t>& B_wt, 
                         word_t D, uint64_t pos_in_B_wt)
        const
        {
            //if (pos_in_B_wt >= B_wt.size()) return;
            if (!(B_wt[pos_in_B_wt] & D)) return;  // aqui debo tener cuidado cuando esté en el ultimo nivel, dado que allí debo acceder al arreglo B original, el del automata. Voy a necesitar una funcion que me de acceso a eso.
            //std::cout << pos_in_B_wt << " " << B_wt[pos_in_B_wt] << " & " << D << std::endl;  
            //if (pos_in_B_wt >= B_wt.size() || !(B_wt[pos_in_B_wt] & D)) return; // no active values in this subtree, get out.

            using std::get;
            if (get<0>(r) > get<1>(r))
                return;
            //is[v.level] = v.offset + get<0>(r);

            if (v.level == m_max_level) {
                res_vec.emplace_back(v.sym);
                cnt_answers += 1; // sdsl::size(r);
                return;
            } /*else {
                rank_off[v.level] = m_tree_rank(is[v.level]);
            }*/
            size_type irb = ilb + (1ULL << (m_max_level-v.level));
            size_type mid = (irb + ilb)>>1;

            auto c_v = expand(v);
            auto c_r = expand(v, r);

            if (!sdsl::empty(get<0>(c_r)) and  mid) {
                _all_active_p_values_in_range(get<0>(c_v),get<0>(c_r),
                                 ilb, is, rank_off,
                                 res_vec, report, cnt_answers, B_wt, D, 2*pos_in_B_wt);
            }
            if (!sdsl::empty(get<1>(c_r))) {
                _all_active_p_values_in_range(get<1>(c_v), get<1>(c_r),
                                 mid, is, rank_off, res_vec, report,
                                 cnt_answers, B_wt, D, 2*pos_in_B_wt+1);
            }
        }

        template<typename word_t>
        void
        all_active_p_values_in_range_test(size_type lb, size_type rb,
                                   std::vector<word_t>& B_wt, word_t D,
                                   std::vector<std::pair<value_type, std::pair<size_type, size_type>>>& res_vec
                                   ) const
        {

            size_type cnt_answers = 0;
            if (lb <= rb) {
                _all_active_p_values_in_range_test<word_t>(root(), {lb, rb}, 0, 
                                     res_vec, /*cnt_answers,*/ B_wt, D, 1);
            }
        }


        template<typename word_t>
        void
        _all_active_p_values_in_range_test(node_type v, range_type r, size_type ilb, 
                         std::vector<std::pair<value_type, std::pair<size_type, size_type>>>& res_vec,
                         /*size_type& cnt_answers,*/ std::vector<word_t>& B_wt,
                         word_t D, uint64_t pos_in_B_wt)
        const
        {
            if (!(B_wt[pos_in_B_wt] & D)) return; 

            using std::get;
            if (get<0>(r) > get<1>(r))
                return;

            if (v.level == m_max_level) {
                res_vec.emplace_back(std::pair<value_type,std::pair<size_type, size_type>>
		                     (
                                      v.sym, std::pair<size_type,size_type>(get<0>(r), get<0>(r)+sdsl::size(r)/*-1*/
				     )));
                //cnt_answers += 1; // sdsl::size(r);
                return;
            }
            
            size_type irb = ilb + (1ULL << (m_max_level-v.level));
            size_type mid = (irb + ilb)>>1;

            auto c_v = expand(v);
            auto c_r = expand(v, r);

            if (!sdsl::empty(get<0>(c_r)) and  mid) {
                _all_active_p_values_in_range_test<word_t>(get<0>(c_v),get<0>(c_r),
                                 ilb,
                                 res_vec, /*report, cnt_answers,*/ B_wt, D, 2*pos_in_B_wt);
            }
            if (!sdsl::empty(get<1>(c_r))) {
                _all_active_p_values_in_range_test<word_t>(get<1>(c_v), get<1>(c_r), 
                                 mid, res_vec, /*report,*/
                                 /*cnt_answers,*/ B_wt, D, 2*pos_in_B_wt+1);
            }
        }


        // Implemented by Diego Arroyuelo
        template <typename word_t>
        std::vector<std::pair<value_type, word_t>>
        all_active_s_values_in_range(size_type lb, size_type rb,
                                   initializable_array<word_t>& D_wt, word_t D,
                                   bool report=true) const
        {

            size_type cnt_answers = 0;
            std::vector<std::pair<value_type, word_t>> res_vec;
            if (lb <= rb) {
                std::vector<size_type> is(m_max_level+1);
                std::vector<size_type> rank_off(m_max_level+1); 
                _all_active_s_values_in_range(root(), {lb, rb}, 0, is,
                                 rank_off, res_vec, report, cnt_answers, D_wt, D, 1);
            }
            return res_vec;
        }

        template<typename word_t>
        //std::vector<std::pair<value_type, std::pair<size_type, size_type>>>
        //std::vector<tuple<value_type, word_t, std::pair<size_type, size_type>>>
        void
        all_active_s_values_in_range_test(size_type lb, size_type rb,
                                   initializable_array<word_t>& D_wt, word_t D,
                                   std::vector<tuple<value_type, word_t, std::pair<size_type, size_type>>>& res_vec,
                                   bool report=true) const
        {
            size_type cnt_answers = 0;
            //std::vector<tuple<value_type, word_t, std::pair<size_type, size_type>>> res_vec;
            if (lb <= rb) {
                //std::vector<size_type> is(m_max_level+1);
                //std::vector<size_type> rank_off(m_max_level+1);
                _all_active_s_values_in_range_test(root(), {lb, rb}, 0, /*is,
                                 rank_off,*/ res_vec, report, cnt_answers, D_wt, D, 1);
            }
            //return res_vec;
        }

        template <typename word_t>
        word_t
        _all_active_s_values_in_range_test(node_type v, range_type r, size_type ilb, /*std::vector<size_type>& is,
                         std::vector<size_type>& rank_off,*/ 
                         std::vector<std::tuple<
                                                value_type,
						word_t,
						std::pair<size_type,size_type>>>& res_vec,
                         bool report, size_type& cnt_answers, initializable_array<word_t>& D_wt,
                         word_t D, uint64_t pos_in_D_wt)
        const
        {
            word_t Dtemp = D_wt.atPos(pos_in_D_wt); 
            if ((Dtemp | D) == Dtemp) return Dtemp;  

            //D_wt[pos_in_D_wt] = /*D_wt.atPos(pos_in_D_wt)*/ Dtemp | (D & ~Dtemp);

            using std::get;
            if (get<0>(r) > get<1>(r))
                return Dtemp;

            if (v.level == m_max_level) {
                D = D & ~Dtemp;//    D_wt.atPos(pos_in_D_wt);
                word_t answ = D_wt[pos_in_D_wt] = /*D_wt.atPos(pos_in_D_wt)*/ Dtemp | D;
                res_vec.emplace_back(std::tuple<value_type, word_t,std::pair<size_type, size_type>>
                                     (v.sym, D,
                                      std::pair<size_type,size_type>(get<0>(r), get<0>(r)+sdsl::size(r))
				     ));
                cnt_answers += 1; // sdsl::size(r);
                return answ;
            } 
            size_type irb = ilb + (1ULL << (m_max_level-v.level));
            size_type mid = (irb + ilb)>>1;

            auto c_v = expand(v);
            auto c_r = expand(v, r);
       
            //D_wt[pos_in_D_wt] = /*D_wt.atPos(pos_in_D_wt)*/ Dtemp | D;
            word_t answ1, answ2;
            if (!sdsl::empty(get<0>(c_r)) and  mid) {
                answ1 = _all_active_s_values_in_range_test(get<0>(c_v),get<0>(c_r),
                                 ilb, res_vec, report, cnt_answers, D_wt, D, 2*pos_in_D_wt);
            } else
                answ1 = Dtemp; //D_wt.atPos(2*pos_in_D_wt);

            if (!sdsl::empty(get<1>(c_r))) {
                answ2 = _all_active_s_values_in_range_test(get<1>(c_v), get<1>(c_r),
                                 mid, res_vec, report,
                                 cnt_answers, D_wt, D, 2*pos_in_D_wt+1);
            } else
                answ2 = Dtemp; // D_wt.atPos(2*pos_in_D_wt+1);

            /*word_t answ =*/ D_wt[pos_in_D_wt] = answ1 & answ2; 
            /*if (D_wt[pos_in_D_wt] != (D_wt.atPos(2*pos_in_D_wt) & D_wt.atPos(2*pos_in_D_wt+1))) {
                cout << "son distintos en posicion " << pos_in_D_wt << endl;
                cout << D_wt[pos_in_D_wt] << " versus " << (D_wt.atPos(2*pos_in_D_wt) & D_wt.atPos(2*pos_in_D_wt+1)) << endl;
                uint64_t j;
                cin >> j;
            }*/
            return answ1 & answ2; //answ;
        }




        // implemented by Diego Arroyuelo
        template <typename word_t>
        void
        _all_active_s_values_in_range(node_type v, range_type r, size_type ilb, std::vector<size_type>& is,
                         std::vector<size_type>& rank_off, std::vector<std::pair<value_type,word_t>>& res_vec,
                         bool report, size_type& cnt_answers, initializable_array<word_t>& D_wt,
                         word_t D, uint64_t pos_in_D_wt)
        const
        {
            //if (pos_in_B_wt >= B_wt.size()) return;
            if ((D_wt.atPos(pos_in_D_wt) | D) == D_wt.atPos(pos_in_D_wt)) return;  // aqui debo tener cuidado cuando esté en el ultimo nivel, dado que allí debo acceder al arreglo B original, el del automata. Voy a necesitar una funcion que me de acceso a eso.
            //std::cout << pos_in_B_wt << " " << B_wt[pos_in_B_wt] << " & " << D << std::endl;
            //if (pos_in_B_wt >= B_wt.size() || !(B_wt[pos_in_B_wt] & D)) return; // no active values in this subtree, get out.

            using std::get;
            if (get<0>(r) > get<1>(r))
                return;
            //is[v.level] = v.offset + get<0>(r);

            if (v.level == m_max_level) {
                D = D & ~D_wt.atPos(pos_in_D_wt);
		D_wt[pos_in_D_wt] = D_wt.atPos(pos_in_D_wt) | D;
                res_vec.emplace_back(std::pair<value_type, word_t>(v.sym,D));
                cnt_answers += 1; // sdsl::size(r);
                return;
            } /*else {
                rank_off[v.level] = m_tree_rank(is[v.level]);
            }*/
            size_type irb = ilb + (1ULL << (m_max_level-v.level));
            size_type mid = (irb + ilb)>>1;

            auto c_v = expand(v);
            auto c_r = expand(v, r);

            if (!sdsl::empty(get<0>(c_r)) and  mid) {
                _all_active_s_values_in_range(get<0>(c_v),get<0>(c_r),
                                 ilb, is, rank_off,
                                 res_vec, report, cnt_answers, D_wt, D, 2*pos_in_D_wt);
            }
            if (!sdsl::empty(get<1>(c_r))) {
                _all_active_s_values_in_range(get<1>(c_v), get<1>(c_r),
                                 mid, is, rank_off, res_vec, report,
                                 cnt_answers, D_wt, D, 2*pos_in_D_wt+1);
            }
            
	    D_wt[pos_in_D_wt] = D_wt.atPos(2*pos_in_D_wt) & D_wt.atPos(2*pos_in_D_wt+1);

        }

        // Implemented by Diego Arroyuelo
        size_type
        rel_min_obj_maj(value_type vlb, value_type vrb, size_type lb) const
        {
            if (vrb > (1ULL << m_max_level))
                vrb = (1ULL << m_max_level);
            if (vlb > vrb)
                return m_size+1;

            size_type ans = m_size+1;
            if (lb < m_size) {
                ans = _rel_min_obj_maj(root(), vlb, vrb, {lb,m_size-1}, 0);
            }
            return ans;
        }

        size_type
        _rel_min_obj_maj(node_type v, value_type vlb, value_type vrb, 
                         range_type r, size_type ilb)
        const
        {
            using std::get;
            //std::cout << "[" << get<0>(r) << "," << get<1>(r) << "]" << std::endl; 
            if (get<0>(r)+1 > get<1>(r)+1)
                return m_size+1;  // OJO con el +1, puede estar de mas

            if (v.level == m_max_level) {
                return get<0>(r); // cuidado con esto
            } 
            size_type irb = ilb + /*std::min(m_sigma-2,*/(size_type)(1ULL << (m_max_level-v.level))/*)*/;
            size_type mid = (irb + ilb)>>1;

            /*std::cout << "[" << ilb << "," << std::min(irb,m_sigma-1) << "] versus [" << vlb << "," << vrb << "]" << std::endl;
            char c;
            std::cin >> c;*/
            
            if (vlb <= ilb and std::min(irb-1,m_sigma-1) <= vrb) {
                /*std::cout << "[" << ilb << "," << std::min(irb-1,m_sigma-1) << "] is included in [" << vlb << "," << vrb << "]" << std::endl;
                char d;
                std::cin >> d;*/
                return get<0>(r);
            }

            std::array<range_type, 2> c_r;
            size_type rnk;
            auto c_v = my_expand(v, r, c_r[0], c_r[1], rnk);

            size_type ans1 = m_size + 1, old_ans1 = m_size + 1, ans2 = m_size + 1;
            if (!sdsl::empty(get<0>(c_r)) and  vlb < mid and mid) {
               old_ans1 =  _rel_min_obj_maj(get<0>(c_v), vlb, std::min(vrb,mid-1), get<0>(c_r), ilb);
               if (old_ans1 != m_size + 1) 
                   ans1 = m_tree_select0(v.offset-rnk + old_ans1 + 1) - v.offset; 
            }
            if (!sdsl::empty(get<1>(c_r)) and vrb >= mid) {
               if (ans1 != m_size + 1) {
                   uint64_t min = (get<1>(get<1>(c_r)) < get<0>(get<1>(c_r))+ans1-old_ans1-1)?get<1>(get<1>(c_r)):get<0>(get<1>(c_r))+ans1-old_ans1-1;
                   ans2 =  _rel_min_obj_maj(get<1>(c_v), std::max(mid, vlb), vrb, {get<0>(get<1>(c_r)), min}, mid);
                   if (ans2 == m_size + 1)
                       return ans1;
                   else 
                       ans2 = m_tree_select1(rnk + ans2 + 1) - v.offset;
               } else {
                   ans2 = _rel_min_obj_maj(get<1>(c_v), std::max(mid, vlb), vrb, get<1>(c_r), mid);
                   if (ans2 != m_size + 1) 
                      ans2 = m_tree_select1(rnk + ans2 + 1) - v.offset;
               }
            }
            
            return (ans1 < ans2)? ans1:ans2;
        }


        size_type
        _rel_min_obj_maj_ant(node_type v, value_type vlb, value_type vrb,
                         range_type r, size_type ilb)
        const
        {
            using std::get;
            /*std::cout << "[" << get<0>(r) << "," << get<1>(r) << "]" << std::endl;
            char a;
            std::cin >> a;*/
            if (get<0>(r)+1 > get<1>(r)+1)
                return m_size+1;  // OJO con el +1, puede estar de mas

            if (v.level == m_max_level) {
                //std::cout << "Retorna " << get<0>(r) << std::endl;
                return get<0>(r); // cuidado con esto
            }
            size_type irb = ilb + (1ULL << (m_max_level-v.level));
            size_type mid = (irb + ilb)>>1;

            auto c_v = expand(v);
            auto c_r = expand(v, r);

            //std::cout << " > Los hijos del nodo actual son [" << get<0>(get<0>(c_r)) << "," << get<1>(get<0>(c_r)) << "] y [" << get<0>(get<1>(c_r)) << "," << get<1>(get<1>(c_r)) << "]" << std::endl;  


            size_type ans1 = m_size + 1, ans2 = m_size + 1;
            size_type rnk = m_size + 1;
            /*std::cout << " > Rango de simbolos = [" << vlb << "," << vrb << "]" << std::endl;
            std::cout << " > Se divide en los rangos [" << vlb << "," << std::min(vrb, mid-1) << "] y [" << std::max(mid,vlb) << "," << vrb << "]" << std::endl;
            std::cout << " > irb = " << irb << std::endl; 
            std::cout << " > mid = " << mid << std::endl;
            */

            if (!sdsl::empty(get<0>(c_r)) and  vlb < mid and mid) {
               //std::cout << "if 1" << std::endl;
               ans1 =  _rel_min_obj_maj_ant(get<0>(c_v), vlb, std::min(vrb,mid-1), get<0>(c_r), ilb);
               //std::cout << ans1 << std::endl;
               if (ans1 != m_size + 1) {
               /*    std::cout << "**************************" << std::endl;
                   std::cout << "ans1 ans1="<< ans1 << std::endl;
                   std::cout << "ans1 offset=" << v.offset << std::endl;
                   std::cout << "ans1 rank=" << m_tree_rank(v.offset) << std::endl;
                   std::cout << "ans1 select=" << m_tree_select0(v.offset-m_tree_rank(v.offset) + ans1 + 1) << std::endl;*/
                   rnk = m_tree_rank(v.offset);
                   ans1 = m_tree_select0(v.offset-rnk/*m_tree_rank(v.offset)*/ + ans1 + 1) - v.offset;
               }
               //std::cout << "ans1=" << ans1 << std::endl;
            } /*else
                std::cout << "No if 1" << std::endl;*/
            if (!sdsl::empty(get<1>(c_r)) and vrb >= mid) {
               //std::cout << "if 2" << std::endl;
               ans2 = _rel_min_obj_maj_ant(get<1>(c_v), std::max(mid, vlb), vrb, get<1>(c_r), mid);
               //std::cout << ans2 << std::endl;
               if (ans2 != m_size + 1) {
                  /* std::cout << "**************************" << std::endl;
                   std::cout << "ans2 ans2="<< ans2 << std::endl;
                   std::cout << "ans2 offset=" << v.offset << std::endl;
                   std::cout << "ans2 rank=" << m_tree_rank(v.offset) << std::endl;
                   std::cout << "ans2 select=" << m_tree_select1(m_tree_rank(v.offset) + ans2 +1) << std::endl;*/
                   if (ans1 == m_size + 1) rnk = m_tree_rank(v.offset);
                   ans2 = m_tree_select1(rnk/*m_tree_rank(v.offset)*/ + ans2 + 1) - v.offset;
               } /*else
                   std::cout << "No if 2" << std::endl;*/
               //std::cout << "ans2=" << ans2 << std::endl;
            }

            return (ans1 < ans2)? ans1:ans2;
        }


        //! range_search_2d searches points in the index interval [lb..rb] and value interval [vlb..vrb].
        /*! \param lb     Left bound of index interval (inclusive)
         *  \param rb     Right bound of index interval (inclusive)
         *  \param vlb    Left bound of value interval (inclusive)
         *  \param vrb    Right bound of value interval (inclusive)
         *  \param report Should the matching points be returned?
         *  \return Pair (#of found points, vector of points), the vector is empty when
         *          report = false.
         */
        std::pair<size_type, std::vector<std::pair<value_type, size_type>>>
        range_search_2d(size_type lb, size_type rb, value_type vlb, value_type vrb,
                        bool report=true) const
        {

            if (vrb > (1ULL << m_max_level))
                vrb = (1ULL << m_max_level);
            if (vlb > vrb)
                return make_pair(0, point_vec_type());
            size_type cnt_answers = 0;
            point_vec_type point_vec;
            if (lb <= rb) {
                std::vector<size_type> is(m_max_level+1);
                std::vector<size_type> rank_off(m_max_level+1);
                _range_search_2d(root(), {lb, rb}, vlb, vrb, 0, is,
                                 rank_off, point_vec, report, cnt_answers);
            }
            return make_pair(cnt_answers, point_vec);
        }

        void
        _range_search_2d(node_type v, range_type r, value_type vlb,
                         value_type vrb, size_type ilb, std::vector<size_type>& is,
                         std::vector<size_type>& rank_off, point_vec_type& point_vec,
                         bool report, size_type& cnt_answers)
        const
        {
            using std::get;
            if (get<0>(r) > get<1>(r))
                return;
            is[v.level] = v.offset + get<0>(r);

            if (v.level == m_max_level) {
                for (size_type j=1; j <= sdsl::size(r) and report; ++j) {
                    size_type i = j;
                    size_type c = v.sym;
                    for (uint32_t k=m_max_level; k>0; --k) {
                        size_type offset = is[k-1];
                        size_type rank_offset = rank_off[k-1];
                        if (c&1) {
                            i = m_tree_select1(rank_offset+i)-offset+1;
                        } else {
                            i = m_tree_select0(offset-rank_offset+i)-offset+1;
                        }
                        c >>= 1;
                    }
                    point_vec.emplace_back(is[0]+i-1, v.sym);
                }
                cnt_answers += sdsl::size(r);
                return;
            } else {
                rank_off[v.level] = m_tree_rank(is[v.level]);
            }
            size_type irb = ilb + (1ULL << (m_max_level-v.level));
            size_type mid = (irb + ilb)>>1;

            auto c_v = expand(v);
            auto c_r = expand(v, r);

            if (!sdsl::empty(get<0>(c_r)) and  vlb < mid and mid) {
                _range_search_2d(get<0>(c_v),get<0>(c_r), vlb,
                                 std::min(vrb,mid-1), ilb, is, rank_off,
                                 point_vec, report, cnt_answers);
            }
            if (!sdsl::empty(get<1>(c_r)) and vrb >= mid) {
                _range_search_2d(get<1>(c_v), get<1>(c_r), std::max(mid, vlb),
                                 vrb, mid, is, rank_off, point_vec, report,
                                 cnt_answers);
            }
        }





    //Implemented by Adrian Gomez-Brandon
    //! count_range_search_2d searches points in the index interval [lb..rb] and value interval [vlb..vrb].
    /*! \param lb     Left bound of index interval (inclusive)
     *  \param rb     Right bound of index interval (inclusive)
     *  \param vlb    Left bound of value interval (inclusive)
     *  \param vrb    Right bound of value interval (inclusive)
     *  \return  #of points in the ranges
     */


    size_type count_range_search_2d(size_type lb, size_type rb,
                                    value_type vlb, value_type vrb) const{

        typedef struct {
            node_type node;
            range_type range;
            size_type ilb;
            size_type irb;
        } nrv_type;

        typedef std::stack<nrv_type> stack_type;


        auto nvrb = (1ULL<< m_max_level)-1;
        vrb = (vrb > nvrb) ? nvrb : vrb;

        std::array<range_type, 2> child_ranges;
        size_type cnt_answers = 0, rnk;
        stack_type stack;
        stack.emplace(nrv_type{root(), {lb, rb}, 0, nvrb});
        while (!stack.empty()) {
            nrv_type x = stack.top(); stack.pop();
            if(is_leaf(x.node) || (vlb <= x.ilb && x.irb <= vrb)) {
                cnt_answers += sdsl::size(x.range);
            }else{

                auto child = my_expand(x.node, x.range, child_ranges[0], child_ranges[1], rnk);
                auto mid   = (x.ilb + x.irb+1)>>1;

                if(!sdsl::empty(get<0>(child_ranges)) && mid && vlb < mid) {
                    stack.emplace(nrv_type{get<0>(child), get<0>(child_ranges), x.ilb, mid - 1});
                }
                if(!sdsl::empty(get<1>(child_ranges)) && vrb >= mid) {
                    stack.emplace(nrv_type{get<1>(child), get<1>(child_ranges), mid, x.irb});
                }
            }
        }
        return cnt_answers;
    }

    size_type count_distinct_values(size_type lb, size_type rb) const{

        typedef struct {
            node_type node;
            range_type range;
        } nr_type;

        typedef std::stack<nr_type> stack_type;

        size_type cnt_answers = 0;
        stack_type stack;
        stack.emplace(nr_type{root(), {lb, rb}});
        while (!stack.empty()) {
            nr_type x = stack.top(); stack.pop();
            if(is_leaf(x.node)) {
                ++cnt_answers;
            }else{
                auto child        = expand(x.node);
                auto child_ranges = expand(x.node, x.range);

                if(!sdsl::empty(get<0>(child_ranges))) {
                    stack.emplace(nr_type{get<0>(child), get<0>(child_ranges)});
                }
                if(!sdsl::empty(get<1>(child_ranges))) {
                    stack.emplace(nr_type{get<1>(child), get<1>(child_ranges)});
                }
            }
        }
        return cnt_answers;
    }

    size_type count_range_search_2d_v2(size_type lb, size_type rb,
                                    value_type vlb, value_type vrb) const{
        if (vrb > (1ULL << m_max_level)) vrb = (1ULL << m_max_level);
        if (vlb > vrb) return 0;
        size_type cnt_answers = 0;
        _count_range_search_2d(root(), {lb, rb}, vlb, vrb, 0,  cnt_answers);
        return cnt_answers;
    }

    void _count_range_search_2d(node_type v, range_type r, const value_type vlb,
                          const value_type vrb, const size_type ilb, size_type& cnt_answers) const {
        using std::get;
        if (get<0>(r) > get<1>(r)) return;
        if (v.level == m_max_level ) { //In a leaf
            cnt_answers += sdsl::size(r);
            return;
        }
        size_type irb = ilb + (1ULL << (m_max_level-v.level));
        //The range of values from the leaves that are descendants of
        // the current node is within [vlb, vrb].
        if(vlb >= ilb && irb <= vrb){
            cnt_answers += sdsl::size(r);
            return;
        }

        //Continue with children nodes
        size_type mid = (irb + ilb)>>1;
        auto c_v = expand(v);
        auto c_r = expand(v, r);

        //Check if the range of the left child is not empty and
        // its range of values intersects with [vlb, vrb]
        if (!sdsl::empty(get<0>(c_r)) && mid && vlb < mid) {
            _count_range_search_2d(get<0>(c_v),get<0>(c_r), vlb, std::min(vrb,mid-1), ilb, cnt_answers);
        }

        //Check if the range of the right child is not empty and
        // its range of values intersects with [vlb, vrb]
        if (!sdsl::empty(get<1>(c_r)) && vrb >= mid) {
            _count_range_search_2d(get<1>(c_v), get<1>(c_r), std::max(mid, vlb), vrb, mid, cnt_answers);
        }
    }

        //! Returns a const_iterator to the first element.
        const_iterator begin()const
        {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const
        {
            return const_iterator(this, size());
        }


        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_sigma, out, child, "sigma");
            written_bytes += m_tree.serialize(out, child, "tree");
            written_bytes += m_tree_rank.serialize(out, child, "tree_rank");
            written_bytes += m_tree_select1.serialize(out, child, "tree_select_1");
            written_bytes += m_tree_select0.serialize(out, child, "tree_select_0");
            written_bytes += write_member(m_max_level, out, child, "max_level");
            written_bytes += m_zero_cnt.serialize(out, child, "zero_cnt");
            written_bytes += m_rank_level.serialize(out, child, "rank_level");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            read_member(m_size, in);
            read_member(m_sigma, in);
            m_tree.load(in);
            m_tree_rank.load(in, &m_tree);
            m_tree_select1.load(in, &m_tree);
            m_tree_select0.load(in, &m_tree);
            read_member(m_max_level, in);
            m_zero_cnt.load(in);
            m_rank_level.load(in);
        }

        //! Represents a node in the wavelet tree
        struct node_type {
            size_type  offset   = 0;
            size_type  size     = 0;
            size_type  level    = 0;
            value_type sym      = 0;

            // Default constructor
            node_type(size_type o=0, size_type sz=0, size_type l=0,
                      value_type sy=0) :
                offset(o), size(sz), level(l), sym(sy) {}

            // Copy constructor
            node_type(const node_type&) = default;

            // Move copy constructor
            node_type(node_type&&) = default;

            // Assignment operator
            node_type& operator=(const node_type&) = default;

            // Move assignment operator
            node_type& operator=(node_type&&) = default;

            // Comparator operator
            bool operator==(const node_type& v) const
            {
                return offset == v.offset;
            }

            // Smaller operator
            bool operator<(const node_type& v) const
            {
                return offset < v.offset;
            }

            // Greater operator
            bool operator>(const node_type& v) const
            {
                return offset > v.offset;
            }
        };


        //! Checks if the node is a leaf node
        bool is_leaf(const node_type& v) const
        {
            return v.level == m_max_level;
        }

        //! Symbol of leaf node v
        value_type sym(const node_type& v) const
        {
            return v.sym;
        }

        //! Random access container to bitvector of node v
        auto bit_vec(const node_type& v) const -> node_bv_container<t_bitvector> {
            return node_bv_container<t_bitvector>(begin(v), end(v));
        }

        //! Random access container to sequence of node v
        auto seq(const node_type& v) const -> random_access_container<std::function<value_type(size_type)>> {
            return random_access_container<std::function<value_type(size_type)>>([&v, this](size_type i)
            {
                node_type vv = v;
                while (!is_leaf(vv)) {
                    auto vs = expand(vv);
                    auto rs = expand(vv, {0, i});
                    bool bit = *(begin(vv)+i);
                    i = std::get<1>(rs[bit]);
                    vv = vs[bit];
                }
                return sym(vv);
            }, size(v));
        }

        //! Indicates if node v is empty
        bool empty(const node_type& v) const
        {
            return v.size == (size_type)0;
        }

        //! Return the size of node v
        auto size(const node_type& v) const -> decltype(v.size)
        {
            return v.size;
        }

        //! Return the root node
        node_type root() const
        {
            return node_type(0, m_size, 0, 0);
        }

        //! Returns the two child nodes of an inner node
        /*! \param v An inner node of a wavelet tree.
         *  \return Return a pair of nodes (left child, right child).
         *  \pre !is_leaf(v)
         */
        std::array<node_type, 2>
        expand(const node_type& v) const
        {
            node_type v_right = v;
            return expand(std::move(v_right));
        }



        std::array<node_type, 2>
        my_expand(const node_type& v, const range_type& r, range_type& left_int, range_type& right_int, size_type& rank_b) const
        {
            node_type v_left, v_right = v;
            rank_b = m_tree_rank(v.offset);
            size_type ones   = m_tree_rank(v.offset+v.size)-rank_b; // ones in [b..size)
            size_type ones_p = rank_b - m_rank_level[v.level];      // ones in [level_b..b)

            auto sp_rank    = m_tree_rank(v.offset + r[0]);
            auto right_size = m_tree_rank(v.offset + r[1] + 1) - sp_rank;
            auto left_size  = (r[1]-r[0]+1)-right_size;
            auto right_sp = sp_rank - rank_b;
            auto left_sp  = r[0] - right_sp;
            
            left_int = {left_sp, left_sp + left_size - 1};
            right_int = {right_sp, right_sp + right_size - 1};

            v_left.offset = (v.level+1)*m_size + (v.offset - v.level*m_size) - ones_p;
            v_left.size   = v.size - ones;
            v_left.level  = v.level + 1;
            v_left.sym    = v.sym<<1;

            v_right.offset = (v.level+1)*m_size + m_zero_cnt[v.level] + ones_p;
            v_right.size   = ones;
            v_right.level  = v.level + 1;
            v_right.sym    = (v.sym<<1)|1;

            return {std::move(v_left), std::move(v_right)};
	}

    //Implemented by Adrian Gomez-Brandon.
    //Given a node v and a set of ranges contained in v, it returns the two child nodes of v and computes the left and
    //right child ranges for each range in the set. If a child range is empty, it is not included in the corresponding vector.
    std::array<node_type, 2>
        my_expand_ranges(const node_type& v, const std::vector<range_type>& ranges,
                         std::vector<range_type>& left_ranges,  std::vector<range_type>& right_ranges) const
        {
            size_type sp_rank, right_size, left_size, right_sp, left_sp;
            node_type v_left, v_right = v;
            size_type rank_b = m_tree_rank(v.offset);
            size_type ones   = m_tree_rank(v.offset+v.size)-rank_b; // ones in [b..size)
            size_type ones_p = rank_b - m_rank_level[v.level];      // ones in [level_b..b)

            for (const auto& r : ranges) {
                sp_rank    = m_tree_rank(v.offset + r[0]);
                right_size = m_tree_rank(v.offset + r[1] + 1) - sp_rank;
                left_size  = (r[1]-r[0]+1)-right_size;
                right_sp = sp_rank - rank_b;
                left_sp  = r[0] - right_sp;

                if (left_size > 0) {
                    left_ranges.push_back({left_sp, left_sp + left_size - 1});
                }
                if (right_size > 0) {
                    right_ranges.push_back({right_sp, right_sp + right_size - 1});
                }

            }

            v_left.offset = (v.level+1)*m_size + (v.offset - v.level*m_size) - ones_p;
            v_left.size   = v.size - ones;
            v_left.level  = v.level + 1;
            v_left.sym    = v.sym<<1;

            v_right.offset = (v.level+1)*m_size + m_zero_cnt[v.level] + ones_p;
            v_right.size   = ones;
            v_right.level  = v.level + 1;
            v_right.sym    = (v.sym<<1)|1;

            return {std::move(v_left), std::move(v_right)};
        }
        

        //! Returns the two child nodes of an inner node
        /*! \param v An inner node of a wavelet tree.
         *  \return Return a pair of nodes (left child, right child).
         *  \pre !is_leaf(v)
         */
        std::array<node_type, 2>
        expand(node_type&& v) const
        {
            node_type v_left;
            size_type rank_b = m_tree_rank(v.offset);
            size_type ones   = m_tree_rank(v.offset+v.size)-rank_b; // ones in [b..size)
            size_type ones_p = rank_b - m_rank_level[v.level];      // ones in [level_b..b)

            v_left.offset = (v.level+1)*m_size + (v.offset - v.level*m_size) - ones_p;
            v_left.size   = v.size - ones;
            v_left.level  = v.level + 1;
            v_left.sym    = v.sym<<1;

            v.offset = (v.level+1)*m_size + m_zero_cnt[v.level] + ones_p;
            v.size   = ones;
            v.level  = v.level + 1;
            v.sym    = (v.sym<<1)|1;

            return {std::move(v_left), v};
        }

        //! Returns for each range its left and right child ranges
        /*! \param v      An inner node of an wavelet tree.
         *  \param ranges A vector of ranges. Each range [s,e]
         *                has to be contained in v=[v_s,v_e].
         *  \return A vector a range pairs. The first element of each
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::array<range_vec_type, 2>
        expand(const node_type& v,
               const range_vec_type& ranges) const
        {
            auto ranges_copy = ranges;
            return expand(v, std::move(ranges_copy));
        }

        //! Returns for each range its left and right child ranges
        /*! \param v      An inner node of an wavelet tree.
         *  \param ranges A vector of ranges. Each range [s,e]
         *                has to be contained in v=[v_s,v_e].
         *  \return A vector a range pairs. The first element of each
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::array<range_vec_type, 2>
        expand(const node_type& v,
               range_vec_type&& ranges) const
        {
            auto v_sp_rank = m_tree_rank(v.offset);  // this is already calculated in expand(v)
            range_vec_type res(ranges.size());
            size_t i = 0;
            for (auto& r : ranges) {
                auto sp_rank    = m_tree_rank(v.offset + r[0]);
                auto right_size = m_tree_rank(v.offset + r[1] + 1)
                                  - sp_rank;
                auto left_size  = (r[1]-r[0]+1)-right_size;

                auto right_sp = sp_rank - v_sp_rank;
                auto left_sp  = r[0] - right_sp;

                r = {left_sp, left_sp + left_size - 1};
                res[i++] = {right_sp, right_sp + right_size - 1};
            }
            return {ranges, std::move(res)};
        }

        //! Returns for a range its left and right child ranges
        /*! \param v An inner node of an wavelet tree.
         *  \param r A ranges [s,e], such that [s,e] is
         *           contained in v=[v_s,v_e].
         *  \return A range pair. The first element of the
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::array<range_type, 2>
        expand(const node_type& v, const range_type& r) const
        {
            auto v_sp_rank = m_tree_rank(v.offset);  // this is already calculated in expand(v)
            auto sp_rank    = m_tree_rank(v.offset + r[0]);
            auto right_size = m_tree_rank(v.offset + r[1] + 1)
                              - sp_rank;
            auto left_size  = (r[1]-r[0]+1)-right_size;

            auto right_sp = sp_rank - v_sp_rank;
            auto left_sp  = r[0] - right_sp;

            return {{{left_sp, left_sp + left_size - 1},
                    {right_sp, right_sp + right_size - 1}
                   }
            };
        }

        //! return the path to the leaf for a given symbol
        std::pair<uint64_t,uint64_t> path(value_type c) const
        {
            return {m_max_level,c};
        }

    private:

        //! Iterator to the begin of the bitvector of inner node v
        auto begin(const node_type& v) const -> decltype(m_tree.begin() + v.offset)
        {
            return m_tree.begin() + v.offset;
        }

        //! Iterator to the begin of the bitvector of inner node v
        auto end(const node_type& v) const -> decltype(m_tree.begin() + v.offset + v.size)
        {
            return m_tree.begin() + v.offset + v.size;
        }
};

}// end namespace sdsl
#endif
