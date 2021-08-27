#include "bit_vector.h"

#include <vector>

#include "gtest/gtest.h"

#include "naive_bit_vector.h"

namespace succinct_bv {

class BitVectorTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    v1_.resize(8, false);
    v1_[0] = true;
    v1_[2] = true;
    v1_[3] = true;
    v1_[7] = true; // 10110001

    v2_.resize(10000, false);
    v2_[111] = true;
    v2_[831] = true;
    v2_[5215] = true;

    v3_.resize(1000000, false);
    n_true_ = 0;

    for (int i = 0; i < v3_.size(); ++i) {
      if (rand() % 2 == 0) {
        ++n_true_;
        v3_[i] = true;
      }
    }

    v4_.resize(1000000, false);
    n_sparse_true_ = 0;

    for (int i = 0; i < v4_.size(); ++i) {
      if (rand() % 1000 == 0) {
        ++n_sparse_true_;
        v4_[i] = true;
      }
    }

    v5_.resize(1000000, false);
    n_mix_true_ = 0;
    bool sparse_mode = true;

    for (int i = 0; i < v5_.size(); ++i) {
      if (i % 10000 == 0) sparse_mode = !sparse_mode;

      if (sparse_mode && rand() % 1000 == 0) {
        ++n_mix_true_;
        v5_[i] = true;
      }

      if (!sparse_mode && rand() % 2 == 0) {
        ++n_mix_true_;;
        v5_[i] = true;
      }
    }
  }

  int n_true_;
  int n_sparse_true_;
  int n_mix_true_;
  std::vector<bool> v1_;
  std::vector<bool> v2_;
  std::vector<bool> v3_;
  std::vector<bool> v4_;
  std::vector<bool> v5_;
  std::deque<bool> d1_ = {false,true,false};
  std::vector<bool> empty = {false,false,false,false,false};
};

TEST_F(BitVectorTest, EmptyWorks) {
    BitVector bv1(empty);
    EXPECT_EQ(32, bv1.Select(0));
    EXPECT_EQ(0, bv1.Rank(0));
}

TEST_F(BitVectorTest, AssignWorks) {
    BitVector bv1(v1_);
    EXPECT_EQ(0, bv1.Select(0));
    EXPECT_EQ(2, bv1.Select(1));
    bv1 = v2_;
    EXPECT_EQ(111u, bv1.Select(0));
    EXPECT_EQ(831u, bv1.Select(1));
    bv1 = d1_;
    EXPECT_EQ(1,bv1.Select(0));

    BitVector bv2(v2_);
    EXPECT_EQ(111u, bv2.Select(0));
    EXPECT_EQ(831u, bv2.Select(1));
    bv2 = v1_;
    EXPECT_EQ(0, bv2.Select(0));
    EXPECT_EQ(2, bv2.Select(1));
    bv2 = d1_;
    EXPECT_EQ(1,bv2.Select(0));

    BitVector bv3;
    try {
        EXPECT_EQ(0, bv3.Select(0));
        EXPECT_EQ(0, bv3.Select(1));
        EXPECT_EQ(0, bv3.Rank(0));
        EXPECT_EQ(0, bv3.Rank(1));
    }
    catch (std::runtime_error) {
        EXPECT_EQ(0,0);
    }
    bv3 = v2_;
    EXPECT_EQ(111u, bv3.Select(0));
    EXPECT_EQ(831u, bv3.Select(1));
}

TEST_F(BitVectorTest, AtWorks) {
    BitVector bv1(v1_);
    EXPECT_EQ(true,bv1.At(0));
    EXPECT_EQ(false,bv1.At(1));
    EXPECT_EQ(true,bv1.At(2));
    EXPECT_EQ(true,bv1.At(3));
    EXPECT_EQ(false,bv1.At(4));
    EXPECT_EQ(false,bv1.At(5));
    EXPECT_EQ(false,bv1.At(6));
    EXPECT_EQ(true,bv1.At(7));

    BitVector bv2(v2_);
    EXPECT_EQ(false, bv2.At(110));
    EXPECT_EQ(true, bv2.At(111));
    EXPECT_EQ(false, bv2.At(112));
    EXPECT_EQ(false, bv2.At(830));
    EXPECT_EQ(true, bv2.At(831));
    EXPECT_EQ(false, bv2.At(832));
    EXPECT_EQ(false, bv2.At(5214));
    EXPECT_EQ(true, bv2.At(5215));
    EXPECT_EQ(false, bv2.At(5216));
    EXPECT_EQ(false, bv2.At(8002));
    EXPECT_EQ(false, bv2.At(304));
    EXPECT_EQ(false, bv2.At(3021));
}

TEST_F(BitVectorTest, RankWorks) {
  BitVector bv1(v1_);

  EXPECT_EQ(1u, bv1.Rank(0));
  EXPECT_EQ(1u, bv1.Rank(1));
  EXPECT_EQ(2u, bv1.Rank(2));
  EXPECT_EQ(3u, bv1.Rank(3));
  EXPECT_EQ(3u, bv1.Rank(4));
  EXPECT_EQ(3u, bv1.Rank(5));
  EXPECT_EQ(3u, bv1.Rank(6));
  EXPECT_EQ(4u, bv1.Rank(7));

  BitVector bv2(v2_);

  EXPECT_EQ(0u, bv2.Rank(110));
  EXPECT_EQ(1u, bv2.Rank(111));
  EXPECT_EQ(1u, bv2.Rank(830));
  EXPECT_EQ(2u, bv2.Rank(831));
  EXPECT_EQ(2u, bv2.Rank(5214));
  EXPECT_EQ(3u, bv2.Rank(5215));
  EXPECT_EQ(3u, bv2.Rank(9999));

  BitVector bv3(v3_);
  NaiveBitVector nbv3(v3_);

  for (int i = 0; i < v3_.size(); ++i)
    EXPECT_EQ(nbv3.Rank(i), bv3.Rank(i));

  BitVector bv4(v4_);
  NaiveBitVector nbv4(v4_);

  for (int i = 0; i < v4_.size(); ++i)
    EXPECT_EQ(nbv4.Rank(i), bv4.Rank(i));

  BitVector bv5(v5_);
  NaiveBitVector nbv5(v5_);

  for (int i = 0; i < v5_.size(); ++i)
    EXPECT_EQ(nbv5.Rank(i), bv5.Rank(i));
}

TEST_F(BitVectorTest, SelectWorks) {
  BitVector bv1(v1_);

  EXPECT_EQ(0, bv1.Select(0));
  EXPECT_EQ(2, bv1.Select(1));
  EXPECT_EQ(3, bv1.Select(2));
  EXPECT_EQ(7, bv1.Select(3));
  //TODO: Need some handling of high indices
  EXPECT_EQ(32, bv1.Select(4));

  BitVector bv2(v2_);

  EXPECT_EQ(111u, bv2.Select(0));
  EXPECT_EQ(831u, bv2.Select(1));
  EXPECT_EQ(5215u, bv2.Select(2));

  BitVector bv3(v3_);
  NaiveBitVector nbv3(v3_);

  for (int i = 0; i < n_true_; ++i)
    EXPECT_EQ(nbv3.Select(i), bv3.Select(i));

  BitVector bv4(v4_);
  NaiveBitVector nbv4(v4_);

  for (int i = 0; i < n_sparse_true_; ++i)
    EXPECT_EQ(nbv4.Select(i), bv4.Select(i));

  BitVector bv5(v5_);
  NaiveBitVector nbv5(v5_);

  for (int i = 0; i < n_mix_true_; ++i)
    EXPECT_EQ(nbv5.Select(i), bv5.Select(i));
}

} // namespace succinct_bv
