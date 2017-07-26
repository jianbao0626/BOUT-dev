/*


 */

#ifndef __SINGLEDATAITERATOR_H__
#define __SINGLEDATAITERATOR_H__

#include <iterator>
#include <iostream>
#include "unused.hxx"
#include <output.hxx>

#ifdef _OPENMP
#include <omp.h>
int SDI_spread_work(int num_work, int thread, int max_thread);
#endif

/*!
 * Set of indices - DataIterator is dereferenced into these
 */
struct SIndices {
  int i;
  int nx, ny, nz;
};

#define DI_GET_END ((void *) NULL)

/*!
 * Provides range-based iteration over indices. 
 * If OpenMP is enabled, then this divides work between threads.
 * 
 * This is used mainly to loop over the indices of fields,
 * and provides convenient ways to index 
 * 
 * Example
 * -------
 * 
 * Start,end values for (x,y,z) can be specified directly:
 * 
 *     for(d = DataIterator(xs, xe, ys, ye, zs, ze); !d.done(); ++d) {
 *       // print index
 *       output.write("%d,%d,%d\n", d.x, d.y, d.z);
 *       // Index into a Field3D variable 'f'
 *       output.write("Value = %e\n", f[d]);
 *     }
 * 
 * Usually DataIterator is used to loop over fields. Field3D::begin()
 * and Field3D::end() return DataIterator objects:
 * 
 *     Field3D f(0.0); // Initialise field
 *     for(auto i : f) { // Loop over all indices, including guard cells
 *       f[i] = i.x; // Indexing using DataIterator
 *     }
 * 
 */
class SingleDataIterator
  : public std::iterator<std::bidirectional_iterator_tag, SIndices> {
private:
  /*!
   * This initialises OpenMP threads if enabled, and
   * divides iteration index ranges between threads
   */
#ifdef _OPENMP
  void omp_init(bool end);
#endif
public:
  /*!
   * Constructor. This sets index ranges.
   * If OpenMP is enabled, the index range is divided
   * between threads using the omp_init method.
   */ 
  SingleDataIterator(int is, int xe, int ye, int ze) : 
#ifndef _OPENMP
    nx(xe), ny(ye), nz(ze),
    i(is), istart(is), imin(istart), iend(nx*ny*nz-1), imax(iend),
#else
    nx(xe), ny(ye), nz(ze),
    imin(is), imax(nx*ny*nz-1), 
#endif
    isEnd(false)
  {
#ifdef _OPENMP
    omp_init(false);
#endif
  }

  /*!
   * set end();
   * use as DataIterator(int,int,int,int,int,int,DI_GET_END);
   */
  SingleDataIterator(int is, int xe, int ye, int ze, void* UNUSED(dummy)) :
#ifndef _OPENMP
    nx(xe), ny(ye), nz(ze),
    i(is), istart(is), imin(istart), iend(nx*ny*nz-1), imax(iend),
#else
    nx(xe), ny(ye), nz(ze),
    imin(is), imax(nx*ny*nz-1), 
#endif
    isEnd(true)
  {
#ifdef _OPENMP
    omp_init(true);
#endif
    next();
  }
  
  /*!
   * The index variables, updated during loop
   * Should make these private and provide getters?
   */
  int i;
  const int nx, ny, nz;

  /// Pre-increment operator. Use this rather than post-increment when possible
  SingleDataIterator& operator++() { next(); return *this; }
  
  /// Post-increment operator
  SingleDataIterator operator++(int) { SingleDataIterator tmp(*this); next(); return tmp; }
  
  /// Pre-decrement operator
  SingleDataIterator& operator--() { prev(); return *this; }
  
  /// Post-decrement operator
  SingleDataIterator operator--(int) { SingleDataIterator tmp(*this); prev(); return tmp; }

  /// Comparison operator. Most common use is in for loops
  inline bool operator!=(const SingleDataIterator& rhs) const {
    //return  !(x == rhs.x && y == rhs.y && z == rhs.z);
    if (rhs.isEnd){
      return !this->done();
    } else {
      return  !(i == rhs.i);
    }
  }
  
  /*!
   * Dereference operators
   * These are needed because the C++11 for loop
   * dereferences the iterator
   */
  SingleDataIterator& operator*() {
    return *this;
  }
  
  /*!
   * Const dereference operator. 
   * Needed because C++11 for loop dereferences the iterator
   */
  const SingleDataIterator& operator*() const {
    return *this;
  }

  /*!
   * Add an offset to the index for general stencils
   */
  const SIndices offset(int dx, int dy, int dz) const {
    if (dz>0){
      int zp=i%nz;
      for (int j=0;j<dz;++j)
        zp=(zp == nz-1 ? 0 : zp+1);
      return { i + ny*nz*dx + nz*dy + zp };
    } else {
      int zm=i%nz;
      for (;dz!= 0;++dz)
        zm = (zm == 0 ? nz-1 : zm-1);
      return { i + ny*nz*dx + nz*dy + zm };
    }
  }
  
  /*
   * Shortcuts for common offsets, one cell
   * in each direction.
   */
  
  /// The index one point +1 in x
  const SIndices xp() const { return { i + ny*nz }; }
  /// The index one point -1 in x
  const SIndices xm() const { return { i - ny*nz }; }
  /// The index one point +1 in y
  const SIndices yp() const { return { i + nz }; }
  /// The index one point -1 in y
  const SIndices ym() const { return { i - nz }; }
  /// The index one point +1 in z. Wraps around zend to zstart
  const SIndices zp() const { return { (i+1)%nz == 0 ? i-nz+1 : i+1}; }
  /// The index one point -1 in z. Wraps around zstart to zend
  const SIndices zm() const { return { i%nz == 0 ? i+nz-1 : i-1 }; }

  /*!
   * Resets DataIterator to the start of the range
   */
  void start() {
    i = istart;
  }

  /*!
   * Resets DataIterator to the start of the range
   */
  void begin() {
    i = istart;
  }

  /*!
   * Sets DataIterator to one index past the end of the range
   */ 
  void end() {
    i = iend;
    next();
  }

  /*!
   * Checks if finished looping. Is this more efficient than
   * using the more idiomatic it != DataIterator::end() ?
   */
  bool done() const {
#ifndef _OPENMP
    return (i > iend) || (i < istart);
#else //_OPENMP
    return (i > iend) || (i < istart);
    //return (x == xend) || x > xend || (x <= xstart)  ;
#endif //_OPENMP
  }
  
private:
  SingleDataIterator(); // Disable null constructor

  //const int nx, ny, nz;
#ifndef _OPENMP
  const int istart;
#else
  int istart;
#endif

  int imin;

#ifndef _OPENMP
  const int iend;
#else
  int iend;
#endif

  int imax;

  const bool isEnd;
  /// Advance to the next index
  void next() {
    ++i;
  }

  /// Rewind to the previous index
  void prev() {
    --i;
  }
};

/*!
 * Specifies a range of indices which can be iterated over
 * and begin() and end() methods for range-based for loops
 * 
 * Example
 * -------
 *
 * Index ranges can be defined manually:
 *
 *     IndexRange r(0, 10, 0, 20, 0, 30);
 *     
 * then iterated over using begin() and end()
 *
 *     for( DataIterator i = r.begin(); i != r.end(); i++ ) {
 *       output.write("%d,%d,%d\n", i.x, i.y, i.z);
 *     }
 *
 * or the more convenient range for loop:
 *
 *     for( auto i : r ) {
 *       output.write("%d,%d,%d\n", i.x, i.y, i.z);
 *     }
 *
 * A common use for this class is to loop over
 * regions of a field:
 *
 *     Field3D f(0.0);
 *     for( auto i : f.region(RGN_NOBNDRY) ) {
 *       f[i] = 1.0;
 *     }
 * 
 * where RGN_NOBNDRY specifies a region not including
 * boundary/guard cells. 
 */
struct SIndexRange {
  int istart, nx, ny, nz;
  int iend = nx*ny*nz-1;
  
  const SingleDataIterator begin() const {
    return SingleDataIterator(istart, nx, ny, nz);
  }
  const SingleDataIterator end() const {
    return SingleDataIterator(istart, nx, ny, nz, DI_GET_END);
  }
};

#ifdef _OPENMP
inline int SDI_spread_work(int work,int cp,int np){
  // Spread work between threads. If number of points do not
  // spread evenly between threads, put the remaining "rest"
  // points on the threads 0 ... rest-1.
  int pp=work/np;
  int rest=work%np;
  int result=pp*cp;
  if (rest > cp){
    result +=cp;
  } else {
    result +=rest;
  }
  return result;
};

inline void SingleDataIterator::omp_init(bool end){
  // In the case of OPENMP we need to calculate the range
  int threads=omp_get_num_threads();
  if (threads > 1){
    int work  = (imax-imin+1);
    int current_thread = omp_get_thread_num();
    int begin = SDI_spread_work(work,current_thread,threads);
    int end   = SDI_spread_work(work,current_thread+1,threads);
    //output << "xmin" << xmin;
    //output << "xmax" << xmax;
    //output << "work" << work;
    //output << "begin" << begin << "\n";
    //output << "end" << end << "\n";
    --end;
    iend   = end;
    istart = begin;
  } else {
    istart = imin;
    iend   = imax;
  }
  if (!end){
    i=istart;
  } else {
    i=iend;
  }
};
#endif

#endif // __SINGLEDATAITERATOR_H__
