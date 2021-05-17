#pragma once
#include <vector>

template <typename T>
class TriangleMatrix {
public:
  TriangleMatrix<T>(size_t n) : n(n), elements((n * (n - 1))/2){}
  void set(size_t i, size_t j, T value) {
    assert(i < n && j < n);
    elements[arrayPos(std::min(i,j), std::max(i, j))] = value;
  }
  T get (size_t i, size_t j){
    assert(i < n && j < n);
    return elements[arrayPos(std::min(i,j), std::max(i, j))];
  }
private:

  //determines position in linear array, make sure that i < j
size_t arrayPos(size_t i, size_t j)  {
    assert(i < j);
    size_t pos = ((i*(2*n-i-1))/2) + (j - i - 1);
    assert(pos < elements.size());
    return pos;
  }

  size_t n;
  std::vector<T> elements;
};
