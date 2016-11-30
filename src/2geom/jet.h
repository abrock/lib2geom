// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: keir@google.com (Keir Mierle)
//
// A simple implementation of N-dimensional dual numbers, for automatically
// computing exact derivatives of functions.
//
// While a complete treatment of the mechanics of automatic differentation is
// beyond the scope of this header (see
// http://en.wikipedia.org/wiki/Automatic_differentiation for details), the
// basic idea is to extend normal arithmetic with an extra element, "e," often
// denoted with the greek symbol epsilon, such that e != 0 but e^2 = 0. Dual
// numbers are extensions of the real numbers analogous to complex numbers:
// whereas complex numbers augment the reals by introducing an imaginary unit i
// such that i^2 = -1, dual numbers introduce an "infinitesimal" unit e such
// that e^2 = 0. Dual numbers have two components: the "real" component and the
// "infinitesimal" component, generally written as x + y*e. Surprisingly, this
// leads to a convenient method for computing exact derivatives without needing
// to manipulate complicated symbolic expressions.
//
// For example, consider the function
//
//   f(x) = x^2 ,
//
// evaluated at 10. Using normal arithmetic, f(10) = 100, and df/dx(10) = 20.
// Next, augument 10 with an infinitesimal to get:
//
//   f(10 + e) = (10 + e)^2
//             = 100 + 2 * 10 * e + e^2
//             = 100 + 20 * e       -+-
//                     --            |
//                     |             +--- This is zero, since e^2 = 0
//                     |
//                     +----------------- This is df/dx!
//
// Note that the derivative of f with respect to x is simply the infinitesimal
// component of the value of f(x + e). So, in order to take the derivative of
// any function, it is only necessary to replace the numeric "object" used in
// the function with one extended with infinitesimals. The class Jet, defined in
// this header, is one such example of this, where substitution is done with
// templates.
//
// To handle derivatives of functions taking multiple arguments, different
// infinitesimals are used, one for each variable to take the derivative of. For
// example, consider a scalar function of two scalar parameters x and y:
//
//   f(x, y) = x^2 + x * y
//
// Following the technique above, to compute the derivatives df/dx and df/dy for
// f(1, 3) involves doing two evaluations of f, the first time replacing x with
// x + e, the second time replacing y with y + e.
//
// For df/dx:
//
//   f(1 + e, y) = (1 + e)^2 + (1 + e) * 3
//               = 1 + 2 * e + 3 + 3 * e
//               = 4 + 5 * e
//
//               --> df/dx = 5
//
// For df/dy:
//
//   f(1, 3 + e) = 1^2 + 1 * (3 + e)
//               = 1 + 3 + e
//               = 4 + e
//
//               --> df/dy = 1
//
// To take the gradient of f with the implementation of dual numbers ("jets") in
// this file, it is necessary to create a single jet type which has components
// for the derivative in x and y, and passing them to a templated version of f:
//
//   template<typename T>
//   T f(const T &x, const T &y) {
//     return x * x + x * y;
//   }
//
//   // The "2" means there should be 2 dual number components.
//   Jet<double, 2> x(0);  // Pick the 0th dual number for x.
//   Jet<double, 2> y(1);  // Pick the 1st dual number for y.
//   Jet<double, 2> z = f(x, y);
//
//   LOG(INFO) << "df/dx = " << z.v[0]
//             << "df/dy = " << z.v[1];
//
// Most users should not use Jet objects directly; a wrapper around Jet objects,
// which makes computing the derivative, gradient, or jacobian of templated
// functors simple, is in autodiff.h. Even autodiff.h should not be used
// directly; instead autodiff_cost_function.h is typically the file of interest.
//
// For the more mathematically inclined, this file implements first-order
// "jets". A 1st order jet is an element of the ring
//
//   T[N] = T[t_1, ..., t_N] / (t_1, ..., t_N)^2
//
// which essentially means that each jet consists of a "scalar" value 'a' from T
// and a 1st order perturbation vector 'v' of length N:
//
//   x = a + \sum_i v[i] t_i
//
// A shorthand is to write an element as x = a + u, where u is the pertubation.
// Then, the main point about the arithmetic of jets is that the product of
// perturbations is zero:
//
//   (a + u) * (b + v) = ab + av + bu + uv
//                     = ab + (av + bu) + 0
//
// which is what operator* implements below. Addition is simpler:
//
//   (a + u) + (b + v) = (a + b) + (u + v).
//
// The only remaining question is how to evaluate the function of a jet, for
// which we use the chain rule:
//
//   f(a + u) = f(a) + f'(a) u
//
// where f'(a) is the (scalar) derivative of f at a.
//
// By pushing these things through sufficiently and suitably templated
// functions, we can do automatic differentiation. Just be sure to turn on
// function inlining and common-subexpression elimination, or it will be very
// slow!
//
// WARNING: Most Ceres users should not directly include this file or know the
// details of how jets work. Instead the suggested method for automatic
// derivatives is to use autodiff_cost_function.h, which is a wrapper around
// both jets.h and autodiff.h to make taking derivatives of cost functions for
// use in Ceres easier.

#ifndef CERES_PUBLIC_JET_H_
#define CERES_PUBLIC_JET_H_

#include <cmath>
#include <iosfwd>
#include <iostream>  // NOLINT
#include <limits>
#include <string>

namespace ceres {

template <typename T = double>
struct Jet {

  // Constructor from scalar: a + 0.
  explicit Jet(const T& value) {
    a = value;
    v = T(0);
  }

  Jet() : a(0), v(0) {}

  // Constructor from scalar plus variable: a + t_i.
  Jet(const T& value, const T& _v) {
    a = value;
    v = _v;
  }

  // Compound operators
  Jet<T>& operator+=(const Jet<T> &y) {
    *this = *this + y;
    return *this;
  }

  Jet<T>& operator-=(const Jet<T> &y) {
    *this = *this - y;
    return *this;
  }

  Jet<T>& operator*=(const Jet<T> &y) {
    *this = *this * y;
    return *this;
  }

  Jet<T>& operator/=(const Jet<T> &y) {
    *this = *this / y;
    return *this;
  }

  // The scalar part.
  T a;

  // The infinitesimal part.
  T v;
};

// Unary +
template<typename T> inline
Jet<T> const& operator+(const Jet<T>& f) {
  return f;
}

// TODO(keir): Try adding __attribute__((always_inline)) to these functions to
// see if it causes a performance increase.

// Unary -
template<typename T> inline
Jet<T> operator-(const Jet<T>&f) {
  return Jet<T>(-f.a, -f.v);
}

// Binary +
template<typename T> inline
Jet<T> operator+(const Jet<T>& f,
                 const Jet<T>& g) {
  return Jet<T>(f.a + g.a, f.v + g.v);
}

// Binary + with a scalar: x + s
template<typename T> inline
Jet<T> operator+(const Jet<T>& f, T s) {
  return Jet<T>(f.a + s, f.v);
}

// Binary + with a scalar: s + x
template<typename T> inline
Jet<T> operator+(T s, const Jet<T>& f) {
  return Jet<T>(f.a + s, f.v);
}

// Binary -
template<typename T> inline
Jet<T> operator-(const Jet<T>& f,
                 const Jet<T>& g) {
  return Jet<T>(f.a - g.a, f.v - g.v);
}

// Binary - with a scalar: x - s
template<typename T> inline
Jet<T> operator-(const Jet<T>& f, T s) {
  return Jet<T>(f.a - s, f.v);
}

// Binary - with a scalar: s - x
template<typename T> inline
Jet<T> operator-(T s, const Jet<T>& f) {
  return Jet<T>(s - f.a, -f.v);
}

// Binary *
template<typename T> inline
Jet<T> operator*(const Jet<T>& f,
                 const Jet<T>& g) {
  return Jet<T>(f.a * g.a, f.a * g.v + f.v * g.a);
}

// Binary * with a scalar: x * s
template<typename T> inline
Jet<T> operator*(const Jet<T>& f, T s) {
  return Jet<T>(f.a * s, f.v * s);
}

// Binary * with a scalar: s * x
template<typename T> inline
Jet<T> operator*(T s, const Jet<T>& f) {
  return Jet<T>(f.a * s, f.v * s);
}

// Binary /
template<typename T> inline
Jet<T> operator/(const Jet<T>& f,
                 const Jet<T>& g) {
  // This uses:
  //
  //   a + u   (a + u)(b - v)   (a + u)(b - v)
  //   ----- = -------------- = --------------
  //   b + v   (b + v)(b - v)        b^2
  //
  // which holds because v*v = 0.
  const T g_a_inverse = T(1.0) / g.a;
  const T f_a_by_g_a = f.a * g_a_inverse;
  return Jet<T>(f.a * g_a_inverse, (f.v - f_a_by_g_a * g.v) * g_a_inverse);
}

// Binary / with a scalar: s / x
template<typename T> inline
Jet<T> operator/(T s, const Jet<T>& g) {
  const T minus_s_g_a_inverse2 = -s / (g.a * g.a);
  return Jet<T>(s / g.a, g.v * minus_s_g_a_inverse2);
}

// Binary / with a scalar: x / s
template<typename T> inline
Jet<T> operator/(const Jet<T>& f, T s) {
  const T s_inverse = T(1.0) / s;
  return Jet<T>(f.a * s_inverse, f.v * s_inverse);
}

// Binary comparison operators for both scalars and jets.
#define CERES_DEFINE_JET_COMPARISON_OPERATOR(op) \
  template<typename T> inline \
  bool operator op(const Jet<T>& f, const Jet<T>& g) { \
  return f.a op g.a; \
} \
  template<typename T> inline \
  bool operator op(const T& s, const Jet<T>& g) { \
  return s op g.a; \
} \
  template<typename T> inline \
  bool operator op(const Jet<T>& f, const T& s) { \
  return f.a op s; \
}
CERES_DEFINE_JET_COMPARISON_OPERATOR( <  )  // NOLINT
CERES_DEFINE_JET_COMPARISON_OPERATOR( <= )  // NOLINT
CERES_DEFINE_JET_COMPARISON_OPERATOR( >  )  // NOLINT
CERES_DEFINE_JET_COMPARISON_OPERATOR( >= )  // NOLINT
CERES_DEFINE_JET_COMPARISON_OPERATOR( == )  // NOLINT
CERES_DEFINE_JET_COMPARISON_OPERATOR( != )  // NOLINT
#undef CERES_DEFINE_JET_COMPARISON_OPERATOR

// Pull some functions from namespace std.
//
// This is necessary because we want to use the same name (e.g. 'sqrt') for
// double-valued and Jet-valued functions, but we are not allowed to put
// Jet-valued functions inside namespace std.
//
// TODO(keir): Switch to "using".
inline double abs     (double x) { return std::abs(x);      }
inline double log     (double x) { return std::log(x);      }
inline double exp     (double x) { return std::exp(x);      }
inline double sqrt    (double x) { return std::sqrt(x);     }
inline double cos     (double x) { return std::cos(x);      }
inline double acos    (double x) { return std::acos(x);     }
inline double sin     (double x) { return std::sin(x);      }
inline double asin    (double x) { return std::asin(x);     }
inline double tan     (double x) { return std::tan(x);      }
inline double atan    (double x) { return std::atan(x);     }
inline double sinh    (double x) { return std::sinh(x);     }
inline double cosh    (double x) { return std::cosh(x);     }
inline double tanh    (double x) { return std::tanh(x);     }
inline double floor   (double x) { return std::floor(x);    }
inline double ceil    (double x) { return std::ceil(x);     }
inline double pow  (double x, double y) { return std::pow(x, y);   }
inline double atan2(double y, double x) { return std::atan2(y, x); }

// In general, f(a + h) ~= f(a) + f'(a) h, via the chain rule.

// abs(x + h) ~= x + h or -(x + h)
template<typename T> inline
Jet<T> abs(const Jet<T>& f) {
  return f.a < T(0.0) ? -f : f;
}

// log(a + h) ~= log(a) + h / a
template<typename T> inline
Jet<T> log(const Jet<T>& f) {
  const T a_inverse = T(1.0) / f.a;
  return Jet<T>(log(f.a), f.v * a_inverse);
}

// exp(a + h) ~= exp(a) + exp(a) h
template<typename T> inline
Jet<T> exp(const Jet<T>& f) {
  const T tmp = exp(f.a);
  return Jet<T>(tmp, tmp * f.v);
}

// sqrt(a + h) ~= sqrt(a) + h / (2 sqrt(a))
template<typename T> inline
Jet<T> sqrt(const Jet<T>& f) {
  const T tmp = sqrt(f.a);
  const T two_a_inverse = T(1.0) / (T(2.0) * tmp);
  return Jet<T>(tmp, f.v * two_a_inverse);
}

// cos(a + h) ~= cos(a) - sin(a) h
template<typename T> inline
Jet<T> cos(const Jet<T>& f) {
  return Jet<T>(cos(f.a), - sin(f.a) * f.v);
}

// acos(a + h) ~= acos(a) - 1 / sqrt(1 - a^2) h
template<typename T> inline
Jet<T> acos(const Jet<T>& f) {
  const T tmp = - T(1.0) / sqrt(T(1.0) - f.a * f.a);
  return Jet<T>(acos(f.a), tmp * f.v);
}

// sin(a + h) ~= sin(a) + cos(a) h
template<typename T> inline
Jet<T> sin(const Jet<T>& f) {
  return Jet<T>(sin(f.a), cos(f.a) * f.v);
}

// asin(a + h) ~= asin(a) + 1 / sqrt(1 - a^2) h
template<typename T> inline
Jet<T> asin(const Jet<T>& f) {
  const T tmp = T(1.0) / sqrt(T(1.0) - f.a * f.a);
  return Jet<T>(asin(f.a), tmp * f.v);
}

// tan(a + h) ~= tan(a) + (1 + tan(a)^2) h
template<typename T> inline
Jet<T> tan(const Jet<T>& f) {
  const T tan_a = tan(f.a);
  const T tmp = T(1.0) + tan_a * tan_a;
  return Jet<T>(tan_a, tmp * f.v);
}

// atan(a + h) ~= atan(a) + 1 / (1 + a^2) h
template<typename T> inline
Jet<T> atan(const Jet<T>& f) {
  const T tmp = T(1.0) / (T(1.0) + f.a * f.a);
  return Jet<T>(atan(f.a), tmp * f.v);
}

// sinh(a + h) ~= sinh(a) + cosh(a) h
template<typename T> inline
Jet<T> sinh(const Jet<T>& f) {
  return Jet<T>(sinh(f.a), cosh(f.a) * f.v);
}

// cosh(a + h) ~= cosh(a) + sinh(a) h
template<typename T> inline
Jet<T> cosh(const Jet<T>& f) {
  return Jet<T>(cosh(f.a), sinh(f.a) * f.v);
}

// tanh(a + h) ~= tanh(a) + (1 - tanh(a)^2) h
template<typename T> inline
Jet<T> tanh(const Jet<T>& f) {
  const T tanh_a = tanh(f.a);
  const T tmp = T(1.0) - tanh_a * tanh_a;
  return Jet<T>(tanh_a, tmp * f.v);
}

// The floor function should be used with extreme care as this operation will
// result in a zero derivative which provides no information to the solver.
//
// floor(a + h) ~= floor(a) + 0
template<typename T> inline
Jet<T> floor(const Jet<T>& f) {
  return Jet<T>(floor(f.a));
}

// The ceil function should be used with extreme care as this operation will
// result in a zero derivative which provides no information to the solver.
//
// ceil(a + h) ~= ceil(a) + 0
template<typename T> inline
Jet<T> ceil(const Jet<T>& f) {
  return Jet<T>(ceil(f.a));
}

// Bessel functions of the first kind with integer order equal to 0, 1, n.
//
// Microsoft has deprecated the j[0,1,n]() POSIX Bessel functions in favour of
// _j[0,1,n]().  Where available on MSVC, use _j[0,1,n]() to avoid deprecated
// function errors in client code (the specific warning is suppressed when
// Ceres itself is built).
inline double BesselJ0(double x) {
#if defined(_MSC_VER) && defined(_j0)
  return _j0(x);
#else
  return j0(x);
#endif
}
inline double BesselJ1(double x) {
#if defined(_MSC_VER) && defined(_j1)
  return _j1(x);
#else
  return j1(x);
#endif
}
inline double BesselJn(int n, double x) {
#if defined(_MSC_VER) && defined(_jn)
  return _jn(n, x);
#else
  return jn(n, x);
#endif
}

// For the formulae of the derivatives of the Bessel functions see the book:
// Olver, Lozier, Boisvert, Clark, NIST Handbook of Mathematical Functions,
// Cambridge University Press 2010.
//
// Formulae are also available at http://dlmf.nist.gov

// See formula http://dlmf.nist.gov/10.6#E3
// j0(a + h) ~= j0(a) - j1(a) h
template<typename T> inline
Jet<T> BesselJ0(const Jet<T>& f) {
  return Jet<T>(BesselJ0(f.a),
                -BesselJ1(f.a) * f.v);
}

// See formula http://dlmf.nist.gov/10.6#E1
// j1(a + h) ~= j1(a) + 0.5 ( j0(a) - j2(a) ) h
template<typename T> inline
Jet<T> BesselJ1(const Jet<T>& f) {
  return Jet<T>(BesselJ1(f.a),
                T(0.5) * (BesselJ0(f.a) - BesselJn(2, f.a)) * f.v);
}

// See formula http://dlmf.nist.gov/10.6#E1
// j_n(a + h) ~= j_n(a) + 0.5 ( j_{n-1}(a) - j_{n+1}(a) ) h
template<typename T> inline
Jet<T> BesselJn(int n, const Jet<T>& f) {
  return Jet<T>(BesselJn(n, f.a),
                T(0.5) * (BesselJn(n - 1, f.a) - BesselJn(n + 1, f.a)) * f.v);
}

// Jet Classification. It is not clear what the appropriate semantics are for
// these classifications. This picks that IsFinite and isnormal are "all"
// operations, i.e. all elements of the jet must be finite for the jet itself
// to be finite (or normal). For IsNaN and IsInfinite, the answer is less
// clear. This takes a "any" approach for IsNaN and IsInfinite such that if any
// part of a jet is nan or inf, then the entire jet is nan or inf. This leads
// to strange situations like a jet can be both IsInfinite and IsNaN, but in
// practice the "any" semantics are the most useful for e.g. checking that
// derivatives are sane.

// The jet is finite if all parts of the jet are finite.
template<typename T> inline
bool IsFinite(const Jet<T>& f) {
  return IsFinite(f.a) && IsFinite(f.v);
}

// The jet is infinite if any part of the jet is infinite.
template<typename T> inline
bool IsInfinite(const Jet<T>& f) {
  return IsInfinite(f.a) || IsInfinite(f.v);
}

// The jet is NaN if any part of the jet is NaN.
template<typename T> inline
bool IsNaN(const Jet<T>& f) {
  return IsNaN(f.a) || IsNaN(f.v);
}

// The jet is normal if all parts of the jet are normal.
template<typename T> inline
bool IsNormal(const Jet<T>& f) {
  return IsNormal(f.a) && IsNormal(f.v);
}

// atan2(b + db, a + da) ~= atan2(b, a) + (- b da + a db) / (a^2 + b^2)
//
// In words: the rate of change of theta is 1/r times the rate of
// change of (x, y) in the positive angular direction.
template<typename T> inline
Jet<T> atan2(const Jet<T>& g, const Jet<T>& f) {
  // Note order of arguments:
  //
  //   f = a + da
  //   g = b + db

  T const tmp = T(1.0) / (f.a * f.a + g.a * g.a);
  return Jet<T>(atan2(g.a, f.a), tmp * (- g.a * f.v + f.a * g.v));
}


// pow -- base is a differentiable function, exponent is a constant.
// (a+da)^p ~= a^p + p*a^(p-1) da
template<typename T> inline
Jet<T> pow(const Jet<T>& f, double g) {
  T const tmp = g * pow(f.a, g - T(1.0));
  return Jet<T>(pow(f.a, g), tmp * f.v);
}

// pow -- base is a constant, exponent is a differentiable function.
// We have various special cases, see the comment for pow(Jet, Jet) for
// analysis:
//
// 1. For f > 0 we have: (f)^(g + dg) ~= f^g + f^g log(f) dg
//
// 2. For f == 0 and g > 0 we have: (f)^(g + dg) ~= f^g
//
// 3. For f < 0 and integer g we have: (f)^(g + dg) ~= f^g but if dg
// != 0, the derivatives are not defined and we return NaN.

template<typename T> inline
Jet<T> pow(double f, const Jet<T>& g) {
  if (f == 0 && g.a > 0) {
    // Handle case 2.
    return Jet<T>(T(0.0));
  }
  if (f < 0 && g.a == floor(g.a)) {
    // Handle case 3.
    Jet<T> ret(pow(f, g.a));
    if (g.v != T(0.0)) {
      // Return a NaN when g.v != 0.
      ret.v = std::numeric_limits<T>::quiet_NaN();
    }
    return ret;
  }
  // Handle case 1.
  T const tmp = pow(f, g.a);
  return Jet<T>(tmp, log(f) * tmp * g.v);
}

// pow -- both base and exponent are differentiable functions. This has a
// variety of special cases that require careful handling.
//
// 1. For f > 0:
//    (f + df)^(g + dg) ~= f^g + f^(g - 1) * (g * df + f * log(f) * dg)
//    The numerical evaluation of f * log(f) for f > 0 is well behaved, even for
//    extremely small values (e.g. 1e-99).
//
// 2. For f == 0 and g > 1: (f + df)^(g + dg) ~= 0
//    This cases is needed because log(0) can not be evaluated in the f > 0
//    expression. However the function f*log(f) is well behaved around f == 0
//    and its limit as f-->0 is zero.
//
// 3. For f == 0 and g == 1: (f + df)^(g + dg) ~= 0 + df
//
// 4. For f == 0 and 0 < g < 1: The value is finite but the derivatives are not.
//
// 5. For f == 0 and g < 0: The value and derivatives of f^g are not finite.
//
// 6. For f == 0 and g == 0: The C standard incorrectly defines 0^0 to be 1
//    "because there are applications that can exploit this definition". We
//    (arbitrarily) decree that derivatives here will be nonfinite, since that
//    is consistent with the behavior for f == 0, g < 0 and 0 < g < 1.
//    Practically any definition could have been justified because mathematical
//    consistency has been lost at this point.
//
// 7. For f < 0, g integer, dg == 0: (f + df)^(g + dg) ~= f^g + g * f^(g - 1) df
//    This is equivalent to the case where f is a differentiable function and g
//    is a constant (to first order).
//
// 8. For f < 0, g integer, dg != 0: The value is finite but the derivatives are
//    not, because any change in the value of g moves us away from the point
//    with a real-valued answer into the region with complex-valued answers.
//
// 9. For f < 0, g noninteger: The value and derivatives of f^g are not finite.

template<typename T> inline
Jet<T> pow(const Jet<T>& f, const Jet<T>& g) {
  if (f.a == 0 && g.a >= 1) {
    // Handle cases 2 and 3.
    if (g.a > 1) {
      return Jet<T>(T(0.0));
    }
    return f;
  }
  if (f.a < 0 && g.a == floor(g.a)) {
    // Handle cases 7 and 8.
    T const tmp = g.a * pow(f.a, g.a - T(1.0));
    Jet<T> ret(pow(f.a, g.a), tmp * f.v);
    if (g.v != T(0.0)) {
      // Return a NaN when g.v != 0.
      ret.v = std::numeric_limits<T>::quiet_NaN();
    }
    return ret;
  }
  // Handle the remaining cases. For cases 4,5,6,9 we allow the log() function
  // to generate -HUGE_VAL or NaN, since those cases result in a nonfinite
  // derivative.
  T const tmp1 = pow(f.a, g.a);
  T const tmp2 = g.a * pow(f.a, g.a - T(1.0));
  T const tmp3 = tmp1 * log(f.a);
  return Jet<T>(tmp1, tmp2 * f.v + tmp3 * g.v);
}

// Note: This has to be in the ceres namespace for argument dependent lookup to
// function correctly. Otherwise statements like CHECK_LE(x, 2.0) fail with
// strange compile errors.
template <typename T>
inline std::ostream &operator<<(std::ostream &s, const Jet<T>& z) {
  s << "[" << z.a << " ; " << z.v << "]";
  return s;
}

}  // namespace ceres

#endif  // CERES_PUBLIC_JET_H_
