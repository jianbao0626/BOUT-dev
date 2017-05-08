///
///

#pragma once

#ifndef __SINGLE_INDEX_ITERATOR_H__
#define __SINGLE_INDEX_ITERATORH__

///
/// @param[in] T is a phantom type parameter. This is so that distinct (strong) types
///            can be created since IterateSingleIndex<int> is a different type to IterateSingleIndex<char>
template<typename T>
struct SingleIndexIterator {
  using type_t = SingleIndexIterator<T>;

  SingleIndexIterator(int ind) : index(ind) {}
  
  int index;

  ///
  /// Dereference operators
  /// These are needed because the C++11 for loop
  /// dereferences the iterator
  ///
  type_t& operator*() {
    return *this;
  }

  ///
  /// Const dereference operator. 
  /// Needed because C++11 for loop dereferences the iterator
  ///
  const type_t& operator*() const {
    return *this;
  }
  
  /// Pre-increment operator
  inline type_t& operator++() {++index; return *this; }

  /// Comparison operator. Most common use is in for loops
  inline bool operator!=(const type_t& rhs) const {
    return index != rhs.index;
  }
};

#endif // __SINGLE_INDEX_ITERATOR_H__

