/*
 * Functions which calculate second derivatives
 */

#include <field3d.hxx>
#include <utils.hxx>

#include "stencil_loops.hxx"

/// Second derivative: Central, 2nd order
struct SECOND_DERIV_C2 {
  inline BoutReal operator()(const Stencil &f) const {
    return f.p + f.m - 2.*f.c;
  }
};

/// Second derivative: Central, 4th order
struct SECOND_DERIV_C4 {
  inline BoutReal operator()(const Stencil &f) const {
    return f.p + f.m - 2.*f.c;
  }
};

/*
 * Second derivatives in X
 */

const Field2D indexD2DX2_C2(const Field2D &f) { 
  loopX1<SECOND_DERIV_C2> op; return op(f); }

const Field3D indexD2DX2_C2(const Field3D &f) { 
  loopX1<SECOND_DERIV_C2> op; return op(f); }

const Field2D indexD2DX2_C4(const Field2D &f) { 
  loopX1<SECOND_DERIV_C2> op; return op(f); }

const Field3D indexD2DX2_C4(const Field3D &f) { 
  loopX1<SECOND_DERIV_C4> op; return op(f); }

/*
 * Second derivatives in Y
 */

const Field2D indexD2DY2_C2(const Field2D &f) { 
  loopY1<SECOND_DERIV_C2> op; return op(f); }

const Field3D indexD2DY2_C2(const Field3D &f) { 
  loopY1<SECOND_DERIV_C2> op; return op(f); }

/*
 * Second derivatives in X
 */

const Field2D indexD2DZ2_C2(const Field2D &f) { 
  loopZ1<SECOND_DERIV_C2> op; return op(f); }

const Field3D indexD2DZ2_C2(const Field3D &f) { 
  loopZ1<SECOND_DERIV_C2> op; return op(f); }

const Field2D indexD2DZ2_C4(const Field2D &f) { 
  loopZ1<SECOND_DERIV_C2> op; return op(f); }

const Field3D indexD2DZ2_C4(const Field3D &f) { 
  loopZ1<SECOND_DERIV_C4> op; return op(f); }

