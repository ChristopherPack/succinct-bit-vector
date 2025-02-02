# Succinct Bit Vector

Succinct bit vector supports following queries in constant time using n + O(nlglgn/lgn) space.

- B: 0,1 vector of length n B[0]B[1]...B[n-1]
- rank(x)
  - number of ones in B[0..x]=B[0]B[1]...B[x].
- select(i)
  - position of i-th 1 from the head.

## Usage
Please build library, include `bit_vector.h` and link library.

`BitVector` can be instantiated from only `deque<bool>` and `vector<bool>`.

`BitVector` can be instantiated with an empty constructor for later assignments.

`bool At(uint64_t x)`, `uint64_t Rank(uint64_t x)` and `uint64_t Select(uint64_t i)` are supported.

`operator=(vector<bool>)` and `operator=(deque<bool>)` are supported.

### Example
```c++
#include <vector>

#include <bit_vector.h>

int main() {
  vector<bool> v(100, false);
  v[10] = true;
  v[50] = true;
  BitVector bv(v);
  // r = 1
  int r = bv.Rank(30);
  // s = 50
  int s = bv.Select(1);
  
  BitVector bv2;
  bv2 = v;
  //r2 = 1
  int r2 = bv2.Rank(30);
}
```

```
$ g++ -std=c++17 -I./succinct_bv/include -L./succinct_bv/lib -lsuccinct_bv example.cc
```

## Build
SSE4.2 and POPCNT is needed.

```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ctest
```

## References
R. Raman, V. Raman, and S. S. Rao. Succinct Indexable Dictionaries with Applications to Encoding k-ary Trees and Multisets, ACM Transactions on Algorithms (TALG) , Vol. 3, Issue 4, 2007.
