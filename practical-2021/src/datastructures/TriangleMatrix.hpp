#pragma once
#include <vector>
#include <iostream>

template <typename T>
class TriangleMatrix {
public:
  TriangleMatrix<T>(size_t n, bool diagonal) : n(n), diagonal(diagonal){
    elements = diagonal ? std::vector<T> ((n * (n + 1))/2) : std::vector<T> ((n * (n - 1))/2);
  }
  void set(size_t i, size_t j, T value) {
    assert((i < n || (i == n && diagonal)) && (j < n || (j == n && diagonal)));
    elements[arrayPos(std::min(i,j), std::max(i, j))] = value;
  }
  T get (size_t i, size_t j){
    assert(i < n && j < n);
    return elements[arrayPos(std::min(i,j), std::max(i, j))];
  }
  //@Softwipe unused function, uncomment for printing
  /*
  void print() {
    for(size_t i = 0; i < n; ++i){
      for(size_t j = i; j < n; ++j){
        if(i == j && ! diagonal) continue;
        std::cout << get(i,j) << "; ";
      }
      std::cout << std::endl;
    }
  }*/
  // @Softwipe unused function
  //std::vector<T> getAsVector() const {return elements;}
private:

  //determines position in linear array, make sure that i < j
size_t arrayPos(size_t i, size_t j)  {
    assert(i < j || (j == i && diagonal));
    size_t pos = diagonal ? (((n*(n+1)) - ((n-i)*(n-i+1)))/2) + (j - i) : ((i*(2*n-i-1))/2) + (j - i - 1);
    assert(pos < elements.size());
    return pos;
  }

  size_t n;
  std::vector<T> elements;
  bool diagonal;
};
