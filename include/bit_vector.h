#ifndef BIT_VECTOR_H_
#define BIT_VECTOR_H_


#include <cstddef>
#include <cstdint>
#include <iostream>

#ifdef _MSC_VER
#define posix_memalign(p, a, s) (((*(p)) = _aligned_malloc((s), (a))), *(p) ?0 :errno)
#endif

#include <deque>
#include <memory>
#include <vector>

namespace succinct_bv {
    class BitVector;
}
void swap(succinct_bv::BitVector& a, succinct_bv::BitVector&);

namespace succinct_bv {
    class BitVector {
    public:
        BitVector() : b_(nullptr) {};

        BitVector(const BitVector& copy);

        BitVector(BitVector&& copy);

        BitVector(const std::deque<bool> &v) : b_(nullptr) { Init(v); }

        BitVector(const std::vector<bool> &v) : b_(nullptr) { Init(v); }

        ~BitVector() {
#ifdef _MSC_VER
            if(this->b_ != nullptr) _aligned_free(b_);
#else
            if(this->b_ != nullptr) free(b_);
#endif
        }

        bool At(uint64_t x) const;

        uint64_t Rank(uint64_t x) const;

        uint64_t Select(uint64_t i) const {
            if (b_ == nullptr) throw std::runtime_error("Bitvector is empty.");
            if(s_.empty()) return 32;
            return s_[i / (64 * 64)]->Select(this, i % (64 * 64));
        }

        size_t n_bytes() const;

        friend void ::swap(BitVector& a, BitVector& b);

        BitVector& operator=(BitVector bv);

        //BitVector& operator=(BitVector& bv);

        BitVector& operator=(BitVector&& copy) noexcept;

        BitVector& operator=(std::vector<bool> const& bv);

        BitVector& operator=(std::vector<bool>&& bv);

        BitVector& operator=(std::deque<bool> const& bv);

        BitVector& operator=(std::deque<bool>&& bv);

    private:

        template<class T> void Init(const T &v) {
            if (v.empty()) {
                throw std::runtime_error("Given container is empty.");
            }

            InitVector(v);
            InitRankIndex();
            InitSelectIndex();
        }

        template<class T> void InitVector(const T &v);

        void InitRankIndex();

        void InitSelectIndex();

        //void InitSelectTable();

        uint8_t SelectOn32bits(uint32_t bits, uint8_t i) const;

        class SelectIndex {
        public:
            virtual ~SelectIndex() = 0;

            virtual uint64_t Select(const BitVector *b, uint16_t i) const = 0;

            virtual size_t n_bytes() const = 0;
        };

        class SelectIndexArray : public SelectIndex {
        public:
            SelectIndexArray(const BitVector *b, const std::vector<uint64_t> &s)
                    : s_(s) {}

            ~SelectIndexArray() override {}

            uint64_t Select(const BitVector *b, uint16_t i) const override {
                return s_[i];
            }

            size_t n_bytes() const override {
                return s_.capacity() * sizeof(uint64_t);
            }

        private:
            std::vector<uint64_t> s_;
        };

        class SelectIndexTree : public SelectIndex {
        public:
            SelectIndexTree(const BitVector *b, const std::vector<uint64_t> &s)
                    : cumsums_(nullptr) { Init(b, s); }

            ~SelectIndexTree() override {
#ifdef _MSC_VER
                if (cumsums_ != nullptr) _aligned_free(cumsums_);
#else
                if (cumsums_ != nullptr) free(cumsums_);
#endif
                 }

            uint64_t Select(const BitVector *b, uint16_t i) const override;

            size_t n_bytes() const override {
                return 8 * n_inner_ * sizeof(uint16_t);
            }

            void Dump() const;

        private:
            void Init(const BitVector *b, const std::vector<uint64_t> &s);

            int height_;
            size_t n_inner_;
            uint64_t first_block_index_;
            uint64_t first_block_offset_;
            /**
             This array stores cumsum of #ones in subtree for each node.
             Cumsum of #ones is at most w^2 and each node has sqrt(w) children.
             However, each node has 8 children in this implementation for efficiency.
             128 bits is needed per node because each node uses 8 int16_t.
             */
            int16_t *cumsums_;
        };

        class SelectTable {
        public:
            SelectTable();
            std::vector<uint8_t> select_table_;
        };

        uint64_t n_b_ = 0;
        // bit vector storing every 1/2 w bits.
        uint32_t *b_;
        // store rank at i * w^2 in the bit vector. rank is at most w bits.
        std::vector<uint64_t> r1_;
        // store rank at i * 1/2 w in a w^2 bits block. rank is at most 2lg(w) bits.
        std::vector<uint16_t> r2_;
        std::vector<std::shared_ptr<SelectIndex> > s_;
        // for compute select, store 8 bits pattern table instead of 1/2 w bits pattern table.
        static inline SelectTable selectTable = SelectTable();
    };
}

#endif // BIT_VECTOR_H_
