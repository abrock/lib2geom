#ifndef BEZIERFIT_H
#define BEZIERFIT_H

#include <2geom/point.h>
#include <2geom/jet.h>
#include <2geom/point.h>
#include <2geom/curve.h>
#include <2geom/bezier.h>
#include <2geom/bezier-curve.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

/**
 * @brief The BezierfitGeomDistanceData struct holds the data necessary for fitting a Bezier curve
 * to
 */
struct BezierfitData {
  /**
     * @brief target Array of target points.
     */
  std::vector<Geom::Point> const& target;
  /**
     * @brief start Start point of the fitted Bezier.
     */
  Geom::Point start;
  /**
     * @brief end End point of the fitted Bezier.
     */
  Geom::Point end;
};

/**
 * @brief The BezierfitGeomDistanceData struct holds the data necessary for fitting a Bezier curve
 * to
 */
struct BezierfitDataFixedTimes {
  /**
     * @brief target Vector of target points.
     */
  std::vector<Geom::Point> const& target;
  /**
   * @brief times Vector of times where the bezier curve should be evaluated
   */
  std::vector<double> const& times;
  /**
     * @brief start Start point of the fitted Bezier.
     */
  Geom::Point start;
  /**
     * @brief end End point of the fitted Bezier.
     */
  Geom::Point end;
};

template<class T>
std::pair<T, T> bezier_templated(
    Geom::Point const start,
    Geom::Point const end,
    T const mid1_x,
    T const mid1_y,
    T const mid2_x,
    T const mid2_y,
    T const t
    ) {

  /* Calculation of bezier(t) - target */
  T const c0 = (T(1) - t) * (T(1) - t) * (T(1) - t);
  T const c1 = T(3) * (T(1) - t) * (T(1) - t) * t;
  T const c2 = T(3) * (T(1) - t) * t * t;
  T const c3 = t * t * t;

  T const bez_x =
      c0 * T(start.x()) +
      c1 * T(mid1_x) +
      c2 * T(mid2_x) +
      c3 * T(end.x());

  T const bez_y =
      c0 * T(start.y()) +
      c1 * T(mid1_y) +
      c2 * T(mid2_y) +
      c3 * T(end.y());

  return std::make_pair(bez_x, bez_y);
}

template<class T>
/**
 * @brief punish_lower_zero
 * @param x
 * @param eps
 * @param scale
 * @return
 */
T punish_lower_zero(T const x, T const eps = T(.1), T const scale = T(1)) {
  if (x >= 0.) {
    return T(0);
  }
  return ceres::sqrt(x * x + eps * eps) * scale;
}

/**
 * @brief BezierfitGeomDistanceF Calculate residuals for GSL nonlinear least squares solver.
 * The target objective is to fit a bezier curve s.t. it has constant distance to a given bezier curve.
 * @param[in] x Input data (mid-point 1, mid-point 2, time points t_i) (adjusted by GSL for each iteration)
 * @param[in] data Constant curve-fitting data: @see BezierfitGeomDistanceData
 * @param[out] f Computed residuals
 * @return
 */
int bezierfit_f (
    const gsl_vector * x,
    void *data,
    gsl_vector * f)
{
  struct BezierfitData * d = static_cast<struct BezierfitData *>(data);
  size_t const n = d->target.size();

  double const mid1_x = gsl_vector_get (x, 0);
  double const mid1_y = gsl_vector_get (x, 1);
  double const mid2_x = gsl_vector_get (x, 2);
  double const mid2_y = gsl_vector_get (x, 3);

  for (size_t ii = 0; ii < n; ii++) {
    double t = gsl_vector_get(x, 4 + ii);
    if (t < 0) {
      gsl_vector_set(f, 2*n, punish_lower_zero(t));
      t = 0;
    }
    else if (t > 1) {
      gsl_vector_set(f, 2*n, punish_lower_zero(1-t));
      t = 1;
    }
    auto const result = bezier_templated(
          d->start,
          d->end,
          mid1_x,
          mid1_y,
          mid2_x,
          mid2_y,
          t
          );
    gsl_vector_set(f, 2*ii, result.first - d->target[ii].x());
    gsl_vector_set(f, 2*ii+1, result.second - d->target[ii].y());
  }
  return GSL_SUCCESS;
}

/**
 * @brief BezierfitGeomDistanceF Calculate residuals for GSL nonlinear least squares solver.
 * The target objective is to fit a bezier curve s.t. it has constant distance to a given bezier curve.
 * @param[in] x Input data (mid-point 1, mid-point 2) (adjusted by GSL for each iteration)
 * @param[in] data Constant curve-fitting data: @see BezierfitGeomDistanceData
 * @param[out] f Computed residuals
 * @return
 */
int bezierfit_fixed_times_f (
    const gsl_vector * x,
    void *data,
    gsl_vector * f)
{
  struct BezierfitDataFixedTimes * d = static_cast<struct BezierfitDataFixedTimes *>(data);
  size_t const n = d->target.size();

  double const mid1_x = gsl_vector_get (x, 0);
  double const mid1_y = gsl_vector_get (x, 1);
  double const mid2_x = gsl_vector_get (x, 2);
  double const mid2_y = gsl_vector_get (x, 3);

  for (size_t ii = 0; ii < n; ii++) {
    double t = d->times[ii];
    if (t < 0) {
      t = 0;
    }
    else if (t > 1) {
      t = 1;
    }
    auto const result = bezier_templated(
          d->start,
          d->end,
          mid1_x,
          mid1_y,
          mid2_x,
          mid2_y,
          t
          );
    gsl_vector_set(f, 2*ii, result.first - d->target[ii].x());
    gsl_vector_set(f, 2*ii+1, result.second - d->target[ii].y());
  }
  return GSL_SUCCESS;
}

/**
 * @brief BezierfitGeomDistanceDF Calculate derivatives of residuals for GSL nonlinear least squares solver.
 * The target objective is to fit a bezier curve s.t. it has constant distance to a given bezier curve.
 * @see BezierfitGeomDistanceF, @see JetDistanceF
 * @param[in] x Input data (mid-point 1, mid-point 2, time points t_i) (adjusted by GSL for each iteration)
 * @param[in] data Constant curve-fitting data: @see BezierfitGeomDistanceData
 * @param[out] f Computed residuals
 * @return
 */
int bezierfit_df (
    const gsl_vector * x,
    void *data,
    gsl_matrix * J)
{
  struct BezierfitData * d = static_cast<struct BezierfitData *>(data);
  size_t const n = d->target.size();

  ceres::Jet<> jets_const[4];
  ceres::Jet<> jets_diff[4];

  for (size_t ii = 0; ii < 4; ++ii) {
    jets_const[ii] = ceres::Jet<>(gsl_vector_get (x, ii), 0.);
    jets_diff[ii]  = ceres::Jet<>(gsl_vector_get (x, ii), 1.);
  }

  gsl_matrix_set_zero(J);
  for (size_t ii = 0; ii < n; ii++) {
    double t = gsl_vector_get(x, 4 + ii);
    if (t < 0) {
      gsl_matrix_set(J, 2*n, 4+ii, punish_lower_zero(ceres::Jet<> (t,1.)).v);
      t = 0;
    }
    else if (t > 1) {
      gsl_matrix_set(J, 2*n, 4+ii, punish_lower_zero(ceres::Jet<> (1.-t,-1.)).v);
      t = 1;
    }
    ceres::Jet<> jet_t(t, 0.);
    for (size_t jj = 0; jj < 4; ++jj) {
      auto const result = bezier_templated(
            d->start,
            d->end,
            (0 == jj) ? jets_diff[0] : jets_const[0],
          (1 == jj) ? jets_diff[1] : jets_const[1],
                                     (2 == jj) ? jets_diff[2] : jets_const[2],
        (3 == jj) ? jets_diff[3] : jets_const[3],
        jet_t
        );
      gsl_matrix_set(J, 2*ii, jj, result.first.v);
      gsl_matrix_set(J, 2*ii+1, jj, result.second.v);
      //std::cout << "(" << result.first.v << ", " << result.second.v << ") ";
    }
    ceres::Jet<> jet_t_diff(t, 1.);
    auto const result = bezier_templated(
          d->start,
          d->end,
          jets_const[0],
        jets_const[1],
        jets_const[2],
        jets_const[3],
        jet_t_diff
        );
    gsl_matrix_set(J, 2*ii, 4+ii, result.first.v);
    gsl_matrix_set(J, 2*ii+1, 4+ii, result.second.v);
    //std::cout << "(" << result.first.v << ", " << result.second.v << ") ";
    //std::cout << std::endl;
  }
  return GSL_SUCCESS;
}

/**
 * @brief BezierfitGeomDistanceDF Calculate derivatives of residuals for GSL nonlinear least squares solver.
 * The target objective is to fit a bezier curve s.t. it has constant distance to a given bezier curve.
 * @see BezierfitGeomDistanceF, @see JetDistanceF
 * @param[in] x Input data (mid-point 1, mid-point 2, time points t_i) (adjusted by GSL for each iteration)
 * @param[in] data Constant curve-fitting data: @see BezierfitGeomDistanceData
 * @param[out] f Computed residuals
 * @return
 */
int bezierfit_fixed_times_df (
    const gsl_vector * x,
    void *data,
    gsl_matrix * J)
{
  struct BezierfitDataFixedTimes * d = static_cast<struct BezierfitDataFixedTimes *>(data);
  size_t const n = d->target.size();

  ceres::Jet<> jets_const[4];
  ceres::Jet<> jets_diff[4];

  for (size_t ii = 0; ii < 4; ++ii) {
    jets_const[ii] = ceres::Jet<>(gsl_vector_get (x, ii), 0.);
    jets_diff[ii]  = ceres::Jet<>(gsl_vector_get (x, ii), 1.);
  }

  gsl_matrix_set_zero(J);
  for (size_t ii = 0; ii < n; ii++) {
    double t = d->times[ii];
    if (t < 0) {
      t = 0;
    }
    else if (t > 1) {
      t = 1;
    }
    ceres::Jet<> jet_t(t, 0.);
    for (size_t jj = 0; jj < 4; ++jj) {
      auto const result = bezier_templated(
            d->start,
            d->end,
            (0 == jj) ? jets_diff[0] : jets_const[0],
          (1 == jj) ? jets_diff[1] : jets_const[1],
                                     (2 == jj) ? jets_diff[2] : jets_const[2],
        (3 == jj) ? jets_diff[3] : jets_const[3],
        jet_t
        );
      gsl_matrix_set(J, 2*ii, jj, result.first.v);
      gsl_matrix_set(J, 2*ii+1, jj, result.second.v);
      //std::cout << "(" << result.first.v << ", " << result.second.v << ") ";
    }
  }
  return GSL_SUCCESS;
}

void
callback_bezier_fit(const size_t iter, void *params,
                    const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double rcond;

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stderr, "iter %2zu: p1_x = %.4f, p1_y = %.4f, p2_x = %.4f, p2_y = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_vector_get(x, 1),
          gsl_vector_get(x, 2),
          gsl_vector_get(x, 3),
          1.0 / rcond,
          gsl_blas_dnrm2(f));

  std::cerr << "t_i: ";
  for (size_t ii = 4; ii < x->size; ++ii) {
    std::cerr << gsl_vector_get(x, ii) << ", ";
  }
  std::cerr << std::endl;


}

/**
 * @brief get_initial_t computes an initial guess for the time points on a bezier curve.
 * @param data[in] Points on the bezier curve for which the times shall be estimated.
 * @param initial_t[out]
 */
std::vector<double> get_initial_t(std::vector<Geom::Point> const& data) {
  std::vector<double> initial_t(data.size(), 0.0);
  if (data.empty()) {
    return initial_t;
  }
  initial_t[0] = 0;
  if (data.size() < 2) {
    return initial_t;
  }
  for (size_t ii = 1; ii < data.size(); ++ii) {
    initial_t[ii] = initial_t[ii-1] + Geom::distance(data[ii-1], data[ii]);
  }
  if (initial_t.back() <= 0) {
    for (size_t ii = 0; ii < data.size(); ++ii) {
      initial_t[ii] = static_cast<double>(ii) / (data.size()-1);
    }
  }
  else {
    for (size_t ii = 0; ii < data.size(); ++ii) {
      initial_t[ii] /= initial_t.back();
    }
  }
  return initial_t;
}

double worst_diff(
    Geom::CubicBezier const& c,
    std::vector<Geom::Point> const& data,
    std::vector<double> const& initial_t
    ) {
  double worst_diff = 0;
  for (size_t ii = 0; ii < data.size(); ++ii) {
    double current_diff = Geom::distance(data[ii], c.pointAt(initial_t[ii]));
    worst_diff = std::max(worst_diff, current_diff);
  }
  return worst_diff;
}

double mean_distance(std::vector<Geom::Point> const& data) {
  if (data.size() < 2) {
    return 1;
  }
  double sum = 0;
  for (size_t ii = 1; ii < data.size(); ++ii) {
    sum += Geom::distance(data[ii-1], data[ii]);
  }
  return sum / (data.size()-1);
}

double left_right_dist(std::vector<Geom::Point> const& data, size_t ii) {
  if (data.empty()) {
    return 0;
  }
  if (data.size() < 2) {
    return 0;
  }
  if (ii+1 >= data.size()) {
    return Geom::distance(data.back(), data[data.size()-2]);
  }
  if (0 == ii) {
    return Geom::distance(data.front(), data[1]);
  }
  return Geom::distance(data[ii-1], data[ii]) + Geom::distance(data[ii], data[ii+1]);
}

void get_initial_guess(Geom::CubicBezier& c, std::vector<Geom::Point> const& data) {
  std::vector<double> initial_t = get_initial_t(data);
  std::string winner = "ig provided by user";
  Geom::CubicBezier copy(c);

  double best_worst_diff = worst_diff(copy, data, initial_t);

  {
    copy.setPoint(1, (2. * copy.initialPoint() + copy.finalPoint())/3.);
    copy.setPoint(2, (copy.initialPoint() + 2. * copy.finalPoint())/3.);
    double const current_worst_diff = worst_diff(copy, data, initial_t);
    if (current_worst_diff < best_worst_diff) {
      best_worst_diff = current_worst_diff;
      c = copy;
      winner = "linear segment";
    }
  }
  if (data.size() > 1) {
    double const dist = distance(c[0], c[3]) / 3.0;
    Geom::Point tang1 = (data[1]             - data.front()).normalized();
    Geom::Point tang2 = (data[data.size()-2] - data.back()).normalized();
    copy.setPoint(1, copy.initialPoint() + dist * tang1);
    copy.setPoint(2, copy.finalPoint() + dist * tang2);
    double const current_worst_diff = worst_diff(copy, data, initial_t);
    if (current_worst_diff < best_worst_diff) {
      best_worst_diff = current_worst_diff;
      c = copy;
      winner = "simple tangents";
    }
  }
  if (data.size() > 1) {
    Geom::Point tang1 = (data[1]             - data.front()) / (initial_t[1]);
    Geom::Point tang2 = (data[data.size()-2] - data.back())  / (1.-initial_t[initial_t.size()-2]);
    copy.setPoint(1, copy.initialPoint() + tang1);
    copy.setPoint(2, copy.finalPoint() + tang2);
    double const current_worst_diff = worst_diff(copy, data, initial_t);
    if (current_worst_diff < best_worst_diff) {
      best_worst_diff = current_worst_diff;
      c = copy;
      winner = "linear tangents";
    }
  }
  if (data.size() > 1) {
    Geom::Point tang1;
    Geom::Point tang2;
    double weight1 = 0, weight2 = 0;
    Geom::Point const front = data.front();
    Geom::Point const back = data.back();
    for (size_t ii = 1; ii+1 < data.size(); ++ii) {
      if (initial_t[ii] > 1e-7) {
        double const t = initial_t[ii];
        double const weight = (1.-t)*(1.-t)*(1.-t);
        tang1 += (data[ii] - front) / t * weight;
        weight1 += weight;
      }
      if (initial_t[ii] < 1-(1e-7)) {
        double const t = initial_t[ii];
        double const weight = t*t*t;
        tang2 += (data[ii] - back) / (1.-t) * weight;
        weight2 += weight;
      }
    }
    copy.setPoint(1, copy.initialPoint() + tang1 / weight1);
    copy.setPoint(2, copy.finalPoint() + tang2 / weight2);
    double const current_worst_diff = worst_diff(copy, data, initial_t);
    if (current_worst_diff < best_worst_diff) {
      best_worst_diff = current_worst_diff;
      c = copy;
      winner = "t^3 tangents";
    }
  }
  for (double scale = .001; scale < 10 && data.size() > 1; scale *= 1.5) {
    Geom::Point tang1;
    Geom::Point tang2;
    double weight1 = 0, weight2 = 0;
    Geom::Point const front = data.front();
    Geom::Point const back = data.back();
    for (size_t ii = 1; ii+1 < data.size(); ++ii) {
      if (initial_t[ii] > 1e-7) {
        double const t = initial_t[ii];
        double const weight = std::exp(-t * scale);
        tang1 += (data[ii] - front) / t * weight;
        weight1 += weight;
      }
      if (initial_t[ii] < 1-(1e-7)) {
        double const t = initial_t[ii];
        double const weight = std::exp(-(1.-t) * scale);
        tang2 += (data[ii] - back) / (1.-t) * weight;
        weight2 += weight;
      }
    }
    copy.setPoint(1, copy.initialPoint() + tang1 / weight1);
    copy.setPoint(2, copy.finalPoint() + tang2 / weight2);
    double const current_worst_diff = worst_diff(copy, data, initial_t);
    if (current_worst_diff < best_worst_diff) {
      best_worst_diff = current_worst_diff;
      c = copy;
      winner = std::string("exp tangents, scale ") + std::to_string(scale);
    }
  }
  for (double scale = .001; scale < 10 && data.size() > 1; scale *= 1.5) {
    Geom::Point tang1;
    Geom::Point tang2;
    double weight1 = 0, weight2 = 0;
    Geom::Point const front = data.front();
    Geom::Point const back = data.back();
    for (size_t ii = 1; ii+1 < data.size(); ++ii) {
      if (initial_t[ii] > 1e-7) {
        double const t = initial_t[ii];
        double const weight = std::exp(- t*t * scale);
        tang1 += (data[ii] - front) / t * weight;
        weight1 += weight;
      }
      if (initial_t[ii] < 1-(1e-7)) {
        double const t = initial_t[ii];
        double const weight = std::exp(- (1.-t)*(1.-t) * scale);
        tang2 += (data[ii] - back) / (1.-t) * weight;
        weight2 += weight;
      }
    }
    copy.setPoint(1, copy.initialPoint() + tang1 / weight1);
    copy.setPoint(2, copy.finalPoint() + tang2 / weight2);
    double const current_worst_diff = worst_diff(copy, data, initial_t);
    if (current_worst_diff < best_worst_diff) {
      best_worst_diff = current_worst_diff;
      c = copy;
      winner = std::string("exp^2 tangents, scale ") + std::to_string(scale);
    }
  }

  for (double scale = .00001; scale < 10 && data.size() > 1; scale *= 1.5) {
    Geom::Point tang1;
    Geom::Point tang2;
    double weight1 = 0, weight2 = 0;
    Geom::Point const front = data.front();
    Geom::Point const back = data.back();
    for (size_t ii = 1; ii+1 < data.size(); ++ii) {
      if (initial_t[ii] > 1e-7) {
        double const t = initial_t[ii];
        double const weight = std::exp(-t * scale) * left_right_dist(data, ii);
        tang1 += (data[ii] - front) / t * weight;
        weight1 += weight;
      }
      if (initial_t[ii] < 1-(1e-7)) {
        double const t = initial_t[ii];
        double const weight = std::exp(-(1.-t) * scale) * left_right_dist(data, ii);
        tang2 += (data[ii] - back) / (1.-t) * weight;
        weight2 += weight;
      }
    }
    copy.setPoint(1, copy.initialPoint() + tang1 / weight1);
    copy.setPoint(2, copy.finalPoint() + tang2 / weight2);
    double const current_worst_diff = worst_diff(copy, data, initial_t);
    if (current_worst_diff < best_worst_diff) {
      best_worst_diff = current_worst_diff;
      c = copy;
      winner = std::string("exp tangents with distance weight, scale ") + std::to_string(scale);
    }
  }
  for (double scale = .00001; scale < 10 && data.size() > 1; scale *= 1.5) {
    Geom::Point tang1;
    Geom::Point tang2;
    double weight1 = 0, weight2 = 0;
    Geom::Point const front = data.front();
    Geom::Point const back = data.back();
    for (size_t ii = 1; ii+1 < data.size(); ++ii) {
      if (initial_t[ii] > 1e-7) {
        double const t = initial_t[ii];
        double const weight = std::exp(- t*t * scale)
            * left_right_dist(data, ii);
        tang1 += (data[ii] - front) / t * weight;
        weight1 += weight;
      }
      if (initial_t[ii] < 1-(1e-7)) {
        double const t = initial_t[ii];
        double const weight = std::exp(- (1.-t)*(1.-t) * scale)
            * left_right_dist(data, ii);
        tang2 += (data[ii] - back) / (1.-t) * weight;
        weight2 += weight;
      }
    }
    copy.setPoint(1, copy.initialPoint() + tang1 / weight1);
    copy.setPoint(2, copy.finalPoint() + tang2 / weight2);
    double const current_worst_diff = worst_diff(copy, data, initial_t);
    if (current_worst_diff < best_worst_diff) {
      best_worst_diff = current_worst_diff;
      c = copy;
      winner = std::string("exp^2 tangents with distance weight, scale ") + std::to_string(scale);
    }
  }

  std::cout << "Best i.g. was: " << winner << std::endl;
}



std::pair<double, double> fit_bezier(Geom::CubicBezier& c, std::vector<Geom::Point> const& data) {
  std::pair<double, double> result(0,0);
  if (data.empty()) {
    return result;
  }
  c.setPoint(0, data.front());
  if (data.size() == 1) {
    c.setPoint(1, data.front());
    c.setPoint(2, data.front());
    c.setPoint(3, data.front());
    return result;
  }
  c.setPoint(3, data.back());
  if (data.size() == 2) {
    c.setPoint(1, (2*data.front() + data.back())/3);
    c.setPoint(2, (data.front() + 2*data.back())/3);
    return result;
  }
  if (data.size() == 3) {
    c.setPoint(1, data[1]);
    c.setPoint(2, data[1]);
    return result;
  }
  if (data.size() == 4) {
    c.setPoint(1, data[1]);
    c.setPoint(2, data[2]);
    return result;
  }
  /*
  c.setPoint(1, (2*data.front() + data.back())/3);
  c.setPoint(2, (data.front() + 2*data.back())/3);
*/
  get_initial_guess(c, data);
  const gsl_multifit_nlinear_type *solver_type = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *workspace;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params =
      gsl_multifit_nlinear_default_parameters();

  /**
     * @brief num_points Number of target points to which the bezier c is fitted
     */
  const size_t num_points = data.size();
  /**
     * @brief num_conditions Number of residuals.
     * For each point we have a x residual and y residual.
     * Additionally we have one residual to discourage values for t outside the range [0,1]
     */
  const size_t num_residuals = 2 * num_points + 1;
  /**
     * @brief num_params Number of free parameters.
     * We have for each target point one parameter t_i (the time on c for which the curve is closest).
     * Furthermore we have the two middle control points.
     */
  const size_t num_params = 4 + num_points;

  struct BezierfitData bezierfit_data = {
    data,
        c.initialPoint(),
        c.finalPoint()
  };
  double param_init[num_params];
  param_init[0] = c.controlPoint(1).x();
  param_init[1] = c.controlPoint(1).y();
  param_init[2] = c.controlPoint(2).x();
  param_init[3] = c.controlPoint(2).y();

  double initial_t[num_points];
  initial_t[0] = 0;
  for (size_t ii = 1; ii < num_points; ++ii) {
    initial_t[ii] = initial_t[ii-1] + Geom::distance(data[ii], data[ii-1]);
  }
  for (size_t ii = 1; ii < num_points; ++ii) {
    initial_t[ii] /= initial_t[num_points-1];
  }
  for (size_t ii = 0; ii < num_points; ++ii) {
    param_init[ii+4] = initial_t[ii];
  }
  gsl_vector_view x = gsl_vector_view_array (param_init, num_params);

  int status, info;

  double const mean_point_dist = mean_distance(data);

  /* Initialize the weights for the fitting */;
  double weights[num_residuals];
  gsl_vector_view wts = gsl_vector_view_array(weights, num_residuals);
  /* First we set every weight to 1 */
  for (size_t ii = 0; ii < num_residuals; ++ii) {
    weights[ii] = 1;
  }
  /* Then we weight individual points by their distance to left and right neighbour
     normalized by the total length of the given piecewise data so this becomes scale-free */
  for (size_t ii = 0; ii < num_points; ++ii) {
    weights[2*ii] = weights[2*ii+1] = left_right_dist(data, ii) / mean_point_dist;
  }



  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 1e-12;

  /* define the function to be minimized */
  fdf.f = bezierfit_f;
  fdf.df = bezierfit_df;   /* set to NULL for finite-difference Jacobian */
  //fdf.df = NULL;   /* set to NULL for finite-difference Jacobian */
  fdf.fvv = NULL;     /* not using geodesic acceleration */
  fdf.n = num_residuals;
  fdf.p = num_params;
  fdf.params = &bezierfit_data;

  /* allocate workspace with default parameters */
  workspace = gsl_multifit_nlinear_alloc (solver_type, &fdf_params, num_residuals, num_params);

  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, workspace);
  //gsl_multifit_nlinear_init (&x.vector, &fdf, workspace);

  /* solve the system with a maximum of 150 iterations */
  status = gsl_multifit_nlinear_driver(150, xtol, gtol, ftol,
                                       callback_bezier_fit, NULL, &info, workspace);

  c.setPoint(1, Geom::Point(gsl_vector_get(workspace->x, 0), gsl_vector_get(workspace->x, 1)));
  c.setPoint(2, Geom::Point(gsl_vector_get(workspace->x, 2), gsl_vector_get(workspace->x, 3)));

  double worst_error = 0;
  double worst_time = 0.5;
  size_t out_of_range = 0;
  for (size_t ii = 0; ii < num_points; ++ii) {
    double const current_error = std::sqrt(
          gsl_vector_get(workspace->f, 2*ii) * gsl_vector_get(workspace->f, 2*ii) +
          gsl_vector_get(workspace->f, 2*ii+1) * gsl_vector_get(workspace->f, 2*ii+1)
          );
    if (worst_error < current_error) {
      worst_error = current_error;
      worst_time = gsl_vector_get(workspace->x, 2 + ii);
    }
    double const current_time = gsl_vector_get(workspace->x, 2 + ii);
    if (current_time < 0 || current_time > 1) {
      out_of_range++;
    }
  }

  gsl_multifit_nlinear_free (workspace);
  //std::cout << out_of_range << " out of range from " << num_points << std::endl;
  return std::make_pair(worst_error, worst_time);
}

std::pair<double, double> fit_bezier_fixed_times(Geom::CubicBezier& c, std::vector<Geom::Point> const& data) {
  std::pair<double, double> result(0,0);
  if (data.empty()) {
    return result;
  }
  c.setPoint(0, data.front());
  if (data.size() == 1) {
    c.setPoint(1, data.front());
    c.setPoint(2, data.front());
    c.setPoint(3, data.front());
    return result;
  }
  c.setPoint(3, data.back());
  if (data.size() == 2) {
    c.setPoint(1, (2*data.front() + data.back())/3);
    c.setPoint(2, (data.front() + 2*data.back())/3);
    return result;
  }
  if (data.size() == 3) {
    c.setPoint(1, data[1]);
    c.setPoint(2, data[1]);
    return result;
  }
  /*
  c.setPoint(1, (2*data.front() + data.back())/3);
  c.setPoint(2, (data.front() + 2*data.back())/3);
*/
  get_initial_guess(c, data);
  const gsl_multifit_nlinear_type *solver_type = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *workspace;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params =
      gsl_multifit_nlinear_default_parameters();

  /**
     * @brief num_points Number of target points to which the bezier c is fitted
     */
  const size_t num_points = data.size();
  /**
     * @brief num_conditions Number of residuals.
     * For each point we have a x residual and y residual.
     */
  const size_t num_residuals = 2 * num_points;
  /**
     * @brief num_params Number of free parameters.
     * We have for each target point one parameter t_i (the time on c for which the curve is closest).
     * Furthermore we have the two middle control points.
     */
  const size_t num_params = 4;

  std::vector<double> initial_t = get_initial_t(data);

  struct BezierfitDataFixedTimes bezierfit_data = {
    data,
        initial_t,
        c.initialPoint(),
        c.finalPoint()
  };
  double param_init[num_params];
  param_init[0] = c.controlPoint(1).x();
  param_init[1] = c.controlPoint(1).y();
  param_init[2] = c.controlPoint(2).x();
  param_init[3] = c.controlPoint(2).y();

  gsl_vector_view x = gsl_vector_view_array (param_init, num_params);

  int status, info;

  double const mean_point_dist = mean_distance(data);

  /* Initialize the weights for the fitting */;
  double weights[num_residuals];
  gsl_vector_view wts = gsl_vector_view_array(weights, num_residuals);
  /* First we set every weight to 1 */
  for (size_t ii = 0; ii < num_residuals; ++ii) {
    weights[ii] = 1;
  }
  /* Then we weight individual points by their distance to left and right neighbour
     normalized by the total length of the given piecewise data so this becomes scale-free */
  for (size_t ii = 0; ii < num_points; ++ii) {
    weights[2*ii] = weights[2*ii+1] = left_right_dist(data, ii) / mean_point_dist;
  }



  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 1e-12;

  /* define the function to be minimized */
  fdf.f = bezierfit_fixed_times_f;
  fdf.df = bezierfit_fixed_times_df;   /* set to NULL for finite-difference Jacobian */
  //fdf.df = NULL;   /* set to NULL for finite-difference Jacobian */
  fdf.fvv = NULL;     /* not using geodesic acceleration */
  fdf.n = num_residuals;
  fdf.p = num_params;
  fdf.params = &bezierfit_data;

  /* allocate workspace with default parameters */
  workspace = gsl_multifit_nlinear_alloc (solver_type, &fdf_params, num_residuals, num_params);

  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, workspace);
  //gsl_multifit_nlinear_init (&x.vector, &fdf, workspace);

  /* solve the system with a maximum of 150 iterations */
  status = gsl_multifit_nlinear_driver(150, xtol, gtol, ftol,
                                       NULL, NULL, &info, workspace);

  c.setPoint(1, Geom::Point(gsl_vector_get(workspace->x, 0), gsl_vector_get(workspace->x, 1)));
  c.setPoint(2, Geom::Point(gsl_vector_get(workspace->x, 2), gsl_vector_get(workspace->x, 3)));

  double worst_error = 0;
  double worst_time = 0.5;
  size_t out_of_range = 0;
  for (size_t ii = 0; ii < num_points; ++ii) {
    double const current_error = std::sqrt(
          gsl_vector_get(workspace->f, 2*ii) * gsl_vector_get(workspace->f, 2*ii) +
          gsl_vector_get(workspace->f, 2*ii+1) * gsl_vector_get(workspace->f, 2*ii+1)
          );
    if (worst_error < current_error) {
      worst_error = current_error;
      worst_time = initial_t[ii];
    }
  }

  gsl_multifit_nlinear_free (workspace);
  //std::cout << out_of_range << " out of range from " << num_points << std::endl;
  return std::make_pair(worst_error, worst_time);
}

#endif // BEZIERFIT_H


