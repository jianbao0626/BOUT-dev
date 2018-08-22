
/*******************************************************************************
 * Limiters
 *******************************************************************************/

/// Van Leer limiter. Used in TVD code
BoutReal VANLEER(BoutReal r) { return r + fabs(r) / (1.0 + fabs(r)); }

// Superbee limiter
BoutReal SUPERBEE(BoutReal r) {
  return BOUTMAX(0.0, BOUTMIN(2. * r, 1.0), BOUTMIN(r, 2.));
}

/*******************************************************************************
 * Basic derivative methods.
 * All expect to have an input grid cell at the same location as the output
 * Hence convert cell centred values -> centred values, or left -> left
 *******************************************************************************/

const BoutReal WENO_SMALL = 1.0e-8; // Small number for WENO schemes

////////////////////// FIRST DERIVATIVES /////////////////////

/// central, 2nd order
BoutReal DDX_C2(stencil &f) { return 0.5 * (f.p - f.m); }

/// central, 4th order
BoutReal DDX_C4(stencil &f) { return (8. * f.p - 8. * f.m + f.mm - f.pp) / 12.; }

/// Central WENO method, 2nd order (reverts to 1st order near shocks)
BoutReal DDX_CWENO2(stencil &f) {
  BoutReal isl, isr, isc;  // Smoothness indicators
  BoutReal al, ar, ac, sa; // Un-normalised weights
  BoutReal dl, dr, dc;     // Derivatives using different stencils

  dc = 0.5 * (f.p - f.m);
  dl = f.c - f.m;
  dr = f.p - f.c;

  isl = SQ(dl);
  isr = SQ(dr);
  isc = (13. / 3.) * SQ(f.p - 2. * f.c + f.m) + 0.25 * SQ(f.p - f.m);

  al = 0.25 / SQ(WENO_SMALL + isl);
  ar = 0.25 / SQ(WENO_SMALL + isr);
  ac = 0.5 / SQ(WENO_SMALL + isc);
  sa = al + ar + ac;

  return (al * dl + ar * dr + ac * dc) / sa;
}

// Smoothing 2nd order derivative
BoutReal DDX_S2(stencil &f) {

  // 4th-order differencing
  BoutReal result = (8. * f.p - 8. * f.m + f.mm - f.pp) / 12.;

  result += SIGN(f.c) * (f.pp - 4. * f.p + 6. * f.c - 4. * f.m + f.mm) / 12.;

  return result;
}

///////////////////// SECOND DERIVATIVES ////////////////////

/// Second derivative: Central, 2nd order
BoutReal D2DX2_C2(stencil &f) { return f.p + f.m - 2. * f.c; }

/// Second derivative: Central, 4th order
BoutReal D2DX2_C4(stencil &f) {
  return (-f.pp + 16. * f.p - 30. * f.c + 16. * f.m - f.mm) / 12.;
}

//////////////////////// UPWIND METHODS ///////////////////////

/// Upwinding: Central, 2nd order
BoutReal VDDX_C2(BoutReal vc, stencil &f) { return vc * 0.5 * (f.p - f.m); }

/// Upwinding: Central, 4th order
BoutReal VDDX_C4(BoutReal vc, stencil &f) {
  return vc * (8. * f.p - 8. * f.m + f.mm - f.pp) / 12.;
}

/// upwind, 1st order
BoutReal VDDX_U1(BoutReal vc, stencil &f) {
  return vc >= 0.0 ? vc * (f.c - f.m) : vc * (f.p - f.c);
}

/// upwind, 2nd order
BoutReal VDDX_U2(BoutReal vc, stencil &f) {
  return vc >= 0.0 ? vc * (1.5 * f.c - 2.0 * f.m + 0.5 * f.mm)
                   : vc * (-0.5 * f.pp + 2.0 * f.p - 1.5 * f.c);
}

/// upwind, 3rd order
BoutReal VDDX_U3(BoutReal vc, stencil &f) {
  return vc >= 0.0 ? vc*(4.*f.p - 12.*f.m + 2.*f.mm + 6.*f.c)/12.
    : vc*(-4.*f.m + 12.*f.p - 2.*f.pp - 6.*f.c)/12.;
}

/// 3rd-order WENO scheme
BoutReal VDDX_WENO3(BoutReal vc, stencil &f) {
  BoutReal deriv, w, r;

  if (vc > 0.0) {
    // Left-biased stencil

    r = (WENO_SMALL + SQ(f.c - 2.0 * f.m + f.mm)) /
        (WENO_SMALL + SQ(f.p - 2.0 * f.c + f.m));
    w = 1.0 / (1.0 + 2.0 * r * r);

    deriv = 0.5 * (f.p - f.m) - 0.5 * w * (-f.mm + 3. * f.m - 3. * f.c + f.p);

  } else {
    // Right-biased

    r = (WENO_SMALL + SQ(f.pp - 2.0 * f.p + f.c)) /
        (WENO_SMALL + SQ(f.p - 2.0 * f.c + f.m));
    w = 1.0 / (1.0 + 2.0 * r * r);

    deriv = 0.5 * (f.p - f.m) - 0.5 * w * (-f.m + 3. * f.c - 3. * f.p + f.pp);
  }

  return vc * deriv;
}

/// 3rd-order CWENO. Uses the upwinding code and split flux
BoutReal DDX_CWENO3(stencil &f) {
  BoutReal a, ma = fabs(f.c);
  // Split flux
  a = fabs(f.m);
  if (a > ma)
    ma = a;
  a = fabs(f.p);
  if (a > ma)
    ma = a;
  a = fabs(f.mm);
  if (a > ma)
    ma = a;
  a = fabs(f.pp);
  if (a > ma)
    ma = a;

  stencil sp, vp, sm, vm;

  sp = f + ma;
  sm = ma - f;

  return VDDX_WENO3(0.5, sp) + VDDX_WENO3(-0.5, sm);
}

//////////////////////// FLUX METHODS ///////////////////////

BoutReal FDDX_U1(stencil &v, stencil &f) {
  // Velocity at lower end
  BoutReal vs = 0.5 * (v.m + v.c);
  BoutReal result = (vs >= 0.0) ? vs * f.m : vs * f.c;
  // and at upper
  vs = 0.5 * (v.c + v.p);
  result -= (vs >= 0.0) ? vs * f.c : vs * f.p;

  return - result;
}

BoutReal FDDX_C2(stencil &v, stencil &f) { return 0.5 * (v.p * f.p - v.m * f.m); }

BoutReal FDDX_C4(stencil &v, stencil &f) {
  return (8. * v.p * f.p - 8. * v.m * f.m + v.mm * f.mm - v.pp * f.pp) / 12.;
}

//////////////////////// MUSCL scheme ///////////////////////

void DDX_KT_LR(const stencil &f, BoutReal &fLp, BoutReal &fRp, BoutReal &fLm,
               BoutReal &fRm) {
  // Limiter functions
  BoutReal phi = SUPERBEE((f.c - f.m) / (f.p - f.c));
  BoutReal phi_m = SUPERBEE((f.m - f.mm) / (f.c - f.m));
  BoutReal phi_p = SUPERBEE((f.p - f.c) / (f.pp - f.p));

  fLp = f.c + 0.5 * phi * (f.p - f.c);
  fRp = f.p - 0.5 * phi_p * (f.pp - f.p);

  fLm = f.m + 0.5 * phi_m * (f.c - f.m);
  fRm = f.c - 0.5 * phi * (f.p - f.c);
}

// du/dt = d/dx(f)  with maximum local velocity Vmax
BoutReal DDX_KT(const stencil &f, const stencil &u, const BoutReal Vmax) {
  BoutReal uLp, uRp, uLm, uRm;
  BoutReal fLp, fRp, fLm, fRm;

  DDX_KT_LR(u, uLp, uRp, uLm, uRm);
  DDX_KT_LR(f, fLp, fRp, fLm, fRm);

  BoutReal Fm = 0.5 * (fRm + fLm - Vmax * (uRm - uLm));
  BoutReal Fp = 0.5 * (fRp + fLp - Vmax * (uRp - uLp));

  return Fm - Fp;
}

/*******************************************************************************
 * Staggered differencing methods
 * These expect the output grid cell to be at a different location to the input
 *
 * The stencil no longer has a value in 'C' (centre)
 * instead, points are shifted as follows:
 *
 * mm  -> -3/2 h
 * m   -> -1/2 h
 * p   -> +1/2 h
 * pp  -? +3/2 h
 *
 * NOTE: Cell widths (dx, dy, dz) are currently defined as centre->centre
 * for the methods above. This is currently not taken account of, so large
 * variations in cell size will cause issues.
 *******************************************************************************/

/////////////////////// FIRST DERIVATIVES //////////////////////
// Map Centre -> Low or Low -> Centre

// Second order differencing (staggered)
BoutReal DDX_C2_stag(stencil &f) { return f.p - f.m; }

BoutReal DDX_C4_stag(stencil &f) { return (27. * (f.p - f.m) - (f.pp - f.mm)) / 24.; }

BoutReal D2DX2_C2_stag(stencil &f) { return (f.pp + f.mm - f.p - f.m) / 2.; }
/////////////////////////// UPWINDING ///////////////////////////
// Map (Low, Centre) -> Centre  or (Centre, Low) -> Low
// Hence v contains only (mm, m, p, pp) fields whilst f has 'c' too
//
// v.p is v at +1/2, v.m is at -1/2

BoutReal VDDX_U1_stag(stencil &v, stencil &f) {
  // Lower cell boundary
  BoutReal result = (v.m >= 0) ? v.m * f.m : v.m * f.c;

  // Upper cell boundary
  result -= (v.p >= 0) ? v.p * f.c : v.p * f.p;

  result *= -1;

  // result is now d/dx(v*f), but want v*d/dx(f) so subtract f*d/dx(v)
  result -= f.c * (v.p - v.m);

  return result;
}

BoutReal VDDX_U2_stag(stencil &v, stencil &f) {
  // Calculate d(v*f)/dx = (v*f)[i+1/2] - (v*f)[i-1/2]

  // Upper cell boundary
  BoutReal result = (v.p >= 0.) ? v.p * (1.5*f.c - 0.5*f.m) : v.p * (1.5*f.p - 0.5*f.pp);

  // Lower cell boundary
  result -= (v.m >= 0.) ? v.m * (1.5*f.m - 0.5*f.mm) : v.m * (1.5*f.c - 0.5*f.p);

  // result is now d/dx(v*f), but want v*d/dx(f) so subtract f*d/dx(v)
  result -= f.c * (v.p - v.m);

  return result;
}

BoutReal VDDX_C2_stag(stencil &v, stencil &f) {
  // Result is needed at location of f: interpolate v to f's location and take an
  // unstaggered derivative of f
  return 0.5 * (v.p + v.m) * 0.5 * (f.p - f.m);
}

BoutReal VDDX_C4_stag(stencil &v, stencil &f) {
  // Result is needed at location of f: interpolate v to f's location and take an
  // unstaggered derivative of f
  return (9. * (v.m + v.p) - v.mm - v.pp) / 16. * (8. * f.p - 8. * f.m + f.mm - f.pp) /
         12.;
}

/////////////////////////// FLUX ///////////////////////////
// Map (Low, Centre) -> Centre  or (Centre, Low) -> Low
// Hence v contains only (mm, m, p, pp) fields whilst f has 'c' too
//
// v.p is v at +1/2, v.m is at -1/2

BoutReal FDDX_U1_stag(stencil &v, stencil &f) {
  // Lower cell boundary
  BoutReal result = (v.m >= 0) ? v.m * f.m : v.m * f.c;

  // Upper cell boundary
  result -= (v.p >= 0) ? v.p * f.c : v.p * f.p;

  return - result;
}



BoutReal D4DX4_C2(stencil &f) { return (f.pp - 4. * f.p + 6. * f.c - 4. * f.m + f.mm); }
