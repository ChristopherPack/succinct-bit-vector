#include "bit_vector.h"

#include <cstdlib>
#include <cmath>

#include <iostream>

//#include <x86intrin.h>
#include <nmmintrin.h>
#include <immintrin.h>

//Uncomment if your processor does not provide
/*uint32_t _tzcnt_u32(uint32_t x) {
    uint32_t tmp = 0;
    uint32_t dst = 0;
    while((tmp < 32) && ((x>>tmp)&1) == 0) {
        tmp++;
        dst++;
    }
}*/

#ifdef _MSC_VER
#define posix_memalign(p, a, s) (((*(p)) = _aligned_malloc((s), (a))), *(p) ?0 :errno)
#endif

using std::vector;
using namespace succinct_bv;

BitVector::BitVector(const BitVector &copy) : b_(nullptr) {
    this->n_b_ = copy.n_b_;
    if(copy.b_ != nullptr) {
        posix_memalign((void **) &b_, 32, n_b_ * sizeof(uint32_t));
        std::copy(copy.b_, copy.b_ + copy.n_b_, this->b_);
    }
    this->r1_.resize(copy.r1_.size());
    std::copy(copy.r1_.begin(),copy.r1_.end(), this->r1_.begin());
    this->r2_.resize(copy.r2_.size());
    std::copy(copy.r2_.begin(),copy.r2_.end(), this->r2_.begin());
    this->s_.resize(copy.s_.size());
    std::copy(copy.s_.begin(),copy.s_.end(), this->s_.begin());
}

BitVector::BitVector(BitVector &&copy) : b_(nullptr) {
    swap(*this, copy);
}

BitVector & BitVector::operator=(BitVector bv) {
    swap(*this, bv);
    return *this;
}

BitVector & BitVector::operator=(BitVector &&copy) noexcept {
    swap(*this, copy);
    return *this;
}

BitVector & BitVector::operator=(std::deque<bool> &&bv) {
#ifdef _MSC_VER
    if(this->b_ != nullptr) _aligned_free(this->b_);
#else
    if(this->b_ != nullptr) free(b_);
#endif
    b_ = nullptr;
    this->n_b_ = 0;
    this->r1_ = {};
    this->r2_ = {};
    this->s_ = {};
    Init(bv);
    return *this;
}

BitVector & BitVector::operator=(std::vector<bool> &&bv) {
#ifdef _MSC_VER
    if(this->b_ != nullptr) _aligned_free(this->b_);
#else
    if(this->b_ != nullptr) free(b_);
#endif
    b_ = nullptr;
    this->n_b_ = 0;
    this->r1_ = {};
    this->r2_ = {};
    this->s_ = {};
    Init(bv);
    return *this;
}

BitVector & BitVector::operator=(const std::deque<bool> &bv) {

#ifdef _MSC_VER
    if(this->b_ != nullptr) _aligned_free(this->b_);
#else
    if(this->b_ != nullptr) free(b_);
#endif
    b_ = nullptr;
    this->n_b_ = 0;
    this->r1_ = {};
    this->r2_ = {};
    this->s_ = {};
    Init(bv);
    return *this;
}

BitVector & BitVector::operator=(const std::vector<bool> &bv) {

#ifdef _MSC_VER
    if(this->b_ != nullptr) _aligned_free(this->b_);
#else
    if(this->b_ != nullptr) free(b_);
#endif
    b_ = nullptr;
    this->n_b_ = 0;
    this->r1_ = {};
    this->r2_ = {};
    this->s_ = {};
    Init(bv);
    return *this;
}

void swap(succinct_bv::BitVector& a, succinct_bv::BitVector& b) {
    using std::swap;
    swap(a.b_,b.b_);
    swap(a.n_b_,b.n_b_);
    swap(a.r1_,b.r1_);
    swap(a.r2_,b.r2_);
    swap(a.s_,b.s_);
}

template<class T>
void BitVector::InitVector(const T &v) {
    uint64_t n = v.size();
    n_b_ = n / 32 + 1;
    posix_memalign((void**)&b_, 32, n_b_ * sizeof(uint32_t));

    if (b_ == nullptr)
        throw std::runtime_error("Could not allocate memory for bit vector.");

    for (uint64_t i = 0; i < n_b_; ++i)
        b_[i] = 0;

    for (uint64_t i = 0; i < n; ++i)
        if (v[i]) b_[i / 32] |= 1u << (32 - 1 - (i % 32));
}

template void BitVector::InitVector<std::deque<bool> >(
const std::deque<bool> &v);

template void BitVector::InitVector<std::vector<bool> >(
const std::vector<bool> &v);

bool BitVector::At(uint64_t x) const {
    if (b_ == nullptr) throw std::runtime_error("Bitvector is empty.");
    return (b_[x / 32] & (1 << (31 - x % 32)));
}

//std::vector<uint8_t> BitVector::select_table_ = std::vector<uint8_t>();

void BitVector::InitRankIndex() {
    // the number of w^2 bits blocks is [n/w^2]+1.
    // every w^2 bits block contains 2w small (1/2 w bits) blocks.
    r1_.reserve(n_b_ / (2 * 64) + 1);
    // the number of 1/2 w bits blocks is [2n/w]+1.
    r2_.reserve(n_b_);

    uint64_t r1_sum = 0;
    uint64_t r2_sum = 0;

    for (uint64_t i = 0; i < n_b_; ++i) {
        if (i % (2 * 64) == 0) {
            r1_.push_back(r1_sum);
            r2_sum = 0;
        }

        r2_.push_back(static_cast<uint16_t>(r2_sum));
        uint64_t count = _mm_popcnt_u32(b_[i]);
        r1_sum += count;
        r2_sum += count;
    }
}

uint64_t BitVector::Rank(uint64_t x) const {
    if (b_ == nullptr) throw std::runtime_error("Bitvector is empty.");
    size_t r2_index = x / 32;
    uint32_t bits = b_[r2_index] >> (32 - 1 - (x % 32));

    // popcnt instruction is used instead of the pattern table for 1/2 w bits.
    return r1_[x / (64 * 64)] + r2_[r2_index] + _mm_popcnt_u32(bits);
}

void BitVector::InitSelectIndex() {
    //InitSelectTable();

    vector<uint64_t> s;
    vector<uint64_t> next_s;

    for (uint64_t i = 0; i < n_b_; ++i) {
        uint8_t count = _mm_popcnt_u32(b_[i]);

        for (uint8_t j = 0; j < count; ++j)
            s.push_back(i * 32 + static_cast<uint64_t>(SelectOn32bits(b_[i], j)));

        if (s.size() > (64 * 64)) {
            for (size_t j = 64 * 64; j < s.size(); ++j)
                next_s.push_back(s[j]);

            s.resize(64 * 64);
        }

        // a block contains w^2 ones.
        if ((s.size() == (64 * 64)) || (i == n_b_ - 1)) {

            // a block is sparse if the size of block > w^4 bits.
            if(!s.empty()) {
                if ((s.back() - s.front() + 1) > 64 * 64 * 64 * 64)
                    s_.push_back(std::make_shared<SelectIndexArray>(this, s));
                else
                    s_.push_back(std::make_shared<SelectIndexTree>(this, s));
            }

            s.clear();
            s.swap(next_s);
        }
    }
}

/*void BitVector::InitSelectTable() {
    if(!select_table_.empty()) {
        return;
    }
    select_table_.resize(8 * 256, 8);

    for (size_t i = 0; i < 256; ++i) {
        uint16_t pattern = static_cast<uint16_t>(i);
        uint16_t index = 0;

        for (uint8_t j = 0; j < 8; ++j) {
            if (pattern & (1u << (8 - 1 - j))) {
                select_table_[(index << 8) + pattern] = j;
                ++index;
            }
        }
    }
}*/

BitVector::SelectTable::SelectTable() {
    if(!select_table_.empty()) {
        return;
    }
    select_table_.resize(8 * 256, 8);

    for (size_t i = 0; i < 256; ++i) {
        uint16_t pattern = static_cast<uint16_t>(i);
        uint16_t index = 0;

        for (uint8_t j = 0; j < 8; ++j) {
            if (pattern & (1u << (8 - 1 - j))) {
                select_table_[(index << 8) + pattern] = j;
                ++index;
            }
        }
    }
}

uint8_t BitVector::SelectOn32bits(uint32_t bits, uint8_t i) const {
    uint8_t j = 0;

    while (j < 4) {
        uint8_t count = _mm_popcnt_u32((bits >> ((3 - j) * 8)) & 0xffU);
        if (i < count) break;
        ++j;
        i -= count;
    }

    return j * 8 + selectTable.select_table_[(i << 8) + ((bits >> ((3 - j) * 8)) & 0xffU)];
}

size_t BitVector::n_bytes() const {
    size_t n = n_b_ * sizeof(uint32_t);
    n += r1_.capacity() * sizeof(uint64_t) + r2_.capacity() * sizeof(uint16_t);

    for (auto &v : s_)
        n += v->n_bytes();

    n += selectTable.select_table_.capacity() * sizeof(uint8_t);

    return n;
}

BitVector::SelectIndex::~SelectIndex() {}

void BitVector::SelectIndexTree::Init(const BitVector *b,
                                      const std::vector<uint64_t> &s) {
    first_block_index_ = s.front() / 32;
    // the first small block in this block may be the last small block in previous block.
    // if so, add offset.
    uint32_t first_block = b->b_[first_block_index_];
    int shift = 32 - s.front() % 32;

    if (shift == 32)
        first_block = 0;
    else
        first_block = b->b_[first_block_index_] >> shift;

    first_block_offset_ = _mm_popcnt_u32(first_block);
    size_t  n_blocks = s.back() / 32 - first_block_index_ + 1;
    size_t  n_generation = 1;
    size_t  n_nodes = 1;
    height_ = 0;

    while (n_generation < n_blocks) {
        ++height_;
        n_generation *= 8;
        n_nodes += n_generation;
    }

    n_inner_ = n_nodes - n_generation;
    std::vector<int16_t> nodes(n_nodes, 0);

    for (size_t i = 0; i < n_blocks; ++i)
        nodes[n_inner_ + i] = _mm_popcnt_u32(b->b_[first_block_index_ + i]);

    size_t start = n_inner_;

    for (int h = height_ - 1; h >= 0; --h) {
        size_t m = static_cast<size_t>(std::pow(8, h));
        start -= m;

        for (size_t i = 0; i < m; ++i) {
            size_t  index = start + i;
            nodes[index] = 0;
            size_t offset = 8 * index + 1;

            for (size_t j = 0; j < 8; ++j) {
                nodes[index] += nodes[offset + j];
            }
        }
    }

    if (n_inner_ == 0) return;

    posix_memalign((void**)&cumsums_, 16, 8 * n_inner_ * sizeof(uint16_t));

    if (cumsums_ == nullptr)
        throw std::runtime_error("Could not allocate memory.");

    unsigned int offset = 0;
    unsigned int count = 1;

    for (int h = 0; h < height_; ++h) {
        for (size_t i = 0; i < count; ++i) {
            size_t node = offset + i;
            size_t offset = 8 * node + 1;
            int16_t sum = 0;

            for (size_t j = 0; j < 8; ++j) {
                sum += nodes[offset + j];
                cumsums_[8 * node + j] = sum;
            }
        }

        offset += count;
        count *= 8;
    }
}

uint64_t BitVector::SelectIndexTree::Select(const BitVector *b, uint16_t i) const
{
    i += static_cast<uint16_t>(first_block_offset_);
    unsigned int node = 0;
    unsigned int child = 0;

    //Following Init -> SelectIndexTree() -> InitSelectIndex, height_ should have an upper bound, since s is at most 64 * 64 bit large. Therefore O(1).
    for (int j = 0; j < height_; ++j) {
        __m128i value = _mm_set_epi16(i, i, i, i, i, i, i, i);
        __m128i *ptr = reinterpret_cast<__m128i*>(&cumsums_[8 * node]);
        __m128i to_child = _mm_load_si128(ptr);
        __m128i cmp = _mm_cmpgt_epi16(to_child, value);
        unsigned int mask = static_cast<unsigned int>(_mm_movemask_epi8(cmp));
        child = _tzcnt_u32(mask) / 2;

        if (child > 0)
            i -= cumsums_[8 * node + child - 1];

        node = 8 * node + child + 1;
    }

    uint64_t block_index = first_block_index_ + node - n_inner_;
    uint64_t l = b->SelectOn32bits(b->b_[block_index], static_cast<uint8_t>(i));

    return l + block_index * 32;
}

void BitVector::SelectIndexTree::Dump() const {
    std::cout << "cumsum" << std::endl;
    size_t offset = 0;
    size_t count = 1;

    for (int h = 0; h < height_; ++h) {
        for (int i = 0; i < count; ++i) {
            size_t node = offset + i;

            for (int j = 0; j < 8; ++j)
                std::cout << cumsums_[8 * node + j] << ", ";

            std::cout << "| ";
        }

        std::cout << std::endl;
        offset += count;
        count *= 8;
    }
}
