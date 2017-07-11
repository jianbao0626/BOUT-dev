/*
 * Testing performance of an iterator over the mesh
 *
 */

#include <bout.hxx>

#include <chrono>
#include <iostream>
#include <iterator>
#include <time.h>
#include "omp.h"

// A simple iterator over a 3D set of indices
class MeshIterator
  : public std::iterator< std::forward_iterator_tag, Indices > {
public:
  /// Constructor. This would set ranges. Could depend on thread number
  MeshIterator() : x(0), y(0), z(0), xstart(0), ystart(0), zstart(0) {
    xend = mesh->LocalNx-1;
    yend = mesh->LocalNy-1;
    zend = mesh->LocalNz;
  }

  MeshIterator(int x, int y, int z) : x(x), y(y), z(z), xstart(0), ystart(0), zstart(0) {
    xend = mesh->LocalNx-1;
    yend = mesh->LocalNy-1;
    zend = mesh->LocalNz;
  }

  /// The index variables, updated during loop
  int x, y, z;

  /// Increment operators
  MeshIterator& operator++() { next(); return *this; }
  MeshIterator& operator++(int) { next(); return *this; }

  // Comparison operator
  bool operator!=(const MeshIterator& rhs) const {
    return (x != rhs.x) || (y != rhs.y) || (z != rhs.z);
  }

  // Dereference operator
  Indices operator*() {
    return {x, y, z};
  }

  /// Checks if finished looping. Is this more efficient than
  /// using the more idiomatic it != MeshIterator::end() ?
  bool isDone() const {
    return x > xend;
  }

private:
  int xstart, xend;
  int ystart, yend;
  int zstart, zend;

  /// Advance to the next index
  void next() {
    z++;
    if(z > zend) {
      z = zstart;
      y++;
      if(y > yend) {
        y = ystart;
        x++;
      }
    }
  }
};

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);

  Field3D a = 1.0;
  Field2D b = 2.0;

  Field3D result;
  result.allocate();

  typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
  typedef std::chrono::duration<double> Duration;
  using namespace std::chrono;
  
  // A single loop over block data
  
  BoutReal *ad = &a(0,0,0);
  BoutReal *bd = &b(0,0);
  BoutReal *rd = &result(0,0,0);
  
  int maxit = 10;

  // force some parallel region
///  int sum;
///  #pragma omp parallel for private(sum) // We don't care correctness of sum 
///                                        // Otherwise, use "reduction(+: sum)"
///  for (int i = 0; i < 100; ++i) {
///    sum += i;
///  }


  // Loop over data so first test doesn't have a disadvantage from caching
  for(int i=0;i<maxit;++i) {
#pragma omp parallel for
    for(int j=0;j<mesh->LocalNx*mesh->LocalNy*mesh->LocalNz;++j) {
      rd[j] = ad[j] + bd[j%mesh->LocalNz];
    }
  }
  
  int len = mesh->LocalNx*mesh->LocalNy*mesh->LocalNz;
  SteadyClock start1 = steady_clock::now();
  for(int i=0;i<maxit;++i) {
#pragma omp parallel for
    for(int j=0;j<len;++j) {
      rd[j] = ad[j] + bd[j%(mesh->LocalNz)];
    }
  }
  Duration elapsed1 = steady_clock::now() - start1;

  SteadyClock start2 = steady_clock::now();
  for(int i=0;i<maxit;++i) {
#pragma omp parallel for
    for(int j=0;j<len;++j) {
      rd[j] = ad[j] + bd[j&(mesh->LocalNz-1)];
    }
  }
  Duration elapsed2 = steady_clock::now() - start2;

  // Nested loops over block data
  SteadyClock start3 = steady_clock::now();
  for(int x=0;x<maxit;x++) {
    for(int i=0;i<mesh->LocalNx;++i) {
      for(int j=0;j<mesh->LocalNy;++j) {
#pragma ivdep
#pragma omp parallel for
        for(int k=0;k<mesh->LocalNz;++k) {
          result(i,j,k) = a(i,j,k) + b(i,j);
        }
      }
    }
  }
  Duration elapsed3 = steady_clock::now() - start3;

  // MeshIterator over block data
  SteadyClock start4 = steady_clock::now();
  for(int x=0;x<maxit;x++) {
#pragma GCC ivdep
    for(MeshIterator i; !i.isDone(); ++i){
      result(i.x,i.y,i.z) = a(i.x,i.y,i.z) + b(i.x,i.y);
    }
  }
  Duration elapsed4 = steady_clock::now() - start4;

  // DataIterator using begin(), end()
  SteadyClock start5 = steady_clock::now();
  for(int x=0;x<10;x++) {
#pragma omp parallel
    {
    for(DataIterator i = result.beginDI(), rend=result.endDI(); i != rend; ++i){
    //for(DataIterator i = result.iterator(); !i.done(); ++i){
      result(i.x,i.y,i.z) = a(i.x,i.y,i.z) + b(i.x,i.y);
    }
    }
}
  Duration elapsed5 = steady_clock::now() - start5;

  // DataIterator with done()
  SteadyClock start6 = steady_clock::now();
  for(int x=0;x<10;x++) {
    //for(DataIterator i = begin(result); !i.done() ; ++i){
#pragma omp parallel
    {
    for(DataIterator i = result.iterator(); !i.done(); ++i){
      result(i.x,i.y,i.z) = a(i.x,i.y,i.z) + b(i.x,i.y);
    }
    }
  }
  Duration elapsed6 = steady_clock::now() - start6;

  // Single index, accessed with &
  SteadyClock start7 = steady_clock::now();
  for(int x=0;x<maxit;x++) {
#pragma GCC ivdep
    for(auto &i : result){
      result[i] = a[i] + b[i];
    }
  }
  Duration elapsed7 = steady_clock::now() - start7;

  // Single index, accessed with %
  SteadyClock start8 = steady_clock::now();
  for (int x=0;x<maxit;++x) {
#pragma omp parallel
{
    for (auto &i : result) {
      result[i] = a[i] + b(i);
    }
  }
}
  Duration elapsed8 = steady_clock::now() - start8;
  
  // Range based DataIterator 
  SteadyClock start9 = steady_clock::now();
  for (int x=0;x<maxit;++x) {
#pragma omp parallel
    {
    for (auto &i : result.region(RGN_ALL)) {
      result[i] = a[i] + b[i];
    }
    }
  }
  Duration elapsed9 = steady_clock::now() - start9;

  SteadyClock start10 = steady_clock::now();
  for (int x=0;x<10;++x) {
#pragma ivdep
#pragma omp parallel
    {
    //for (const auto &i : result) {
    for(SingleDataIterator i = result.Siterator(); !i.done(); ++i){
      result(i) = a(i) + b(i);
    }
    }
  }
  Duration elapsed10 = steady_clock::now() - start10;
  
  output << "TIMING\n======\n";
  output << "C loop %                   : " << elapsed1.count() << std::endl;
  output << "C loop &                   : " << elapsed2.count() << std::endl;
  output << "----- (x,y,z) indexing ----" << std::endl;
  output << "Nested loops               : " << elapsed3.count() << std::endl;
  output << "MeshIterator               : " << elapsed4.count() << std::endl;
  output << "DataIterator (begin/end)   : " << elapsed5.count() << std::endl;
  output << "DataIterator (begin/done)  : " << elapsed6.count() << std::endl;
  output << "------ [i] indexing -------" << std::endl;
  output << "Single index &             : " << elapsed7.count() << std::endl;
  output << "Single index %             : " << elapsed8.count() << std::endl;
  output << "Three indices              : " << elapsed9.count() << std::endl;
  output << "SingleDataIterator         : " << elapsed10.count() << std::endl;

  BoutFinalise();
  return 0;
}
