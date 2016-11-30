/** @file
 * @brief Unit tests for Affine
 * Uses the Google Testing Framework
 *//*
 * Authors:
 *   Alexander Brock <>
 *
 * Copyright 2010 Authors
 *
 * This library is free software; you can redistribute it and/or
 * modify it either under the terms of the GNU Lesser General Public
 * License version 2.1 as published by the Free Software Foundation
 * (the "LGPL") or, at your option, under the terms of the Mozilla
 * Public License Version 1.1 (the "MPL"). If you do not alter this
 * notice, a recipient may use your version of this file under either
 * the MPL or the LGPL.
 *
 * You should have received a copy of the LGPL along with this library
 * in the file COPYING-LGPL-2.1; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * You should have received a copy of the MPL along with this library
 * in the file COPYING-MPL-1.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * This software is distributed on an "AS IS" basis, WITHOUT WARRANTY
 * OF ANY KIND, either express or implied. See the LGPL or the MPL for
 * the specific language governing rights and limitations.
 */
#include <chrono>

#include <gtest/gtest.h>
#include <2geom/affine.h>
#include <2geom/transforms.h>
#include <2geom/svg-path-parser.h>
#include <2geom/svg-path-writer.h>
#include <2geom/sbasis-to-bezier.h>

#include "libRunningStats/runningstats.h"
#include <boost/filesystem.hpp>
#include "../toys/helper/geom-pathstroke.h"

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

/* number of data points to fit */
#define N 40

struct data {
    size_t n;
    double * y;
};



int
expb_f (const gsl_vector * x, void *data,
        gsl_vector * f)
{
    size_t n = ((struct data *)data)->n;
    double *y = ((struct data *)data)->y;

    double A = gsl_vector_get (x, 0);
    double lambda = gsl_vector_get (x, 1);
    double b = gsl_vector_get (x, 2);

    size_t i;

    for (i = 0; i < n; i++)
    {
        /* Model Yi = A * exp(-lambda * i) + b */
        double t = i;
        double Yi = A * exp (-lambda * t) + b;
        gsl_vector_set (f, i, Yi - y[i]);
    }

    return GSL_SUCCESS;
}


int
expb_df (const gsl_vector * x, void *data,
         gsl_matrix * J)
{
    size_t n = ((struct data *)data)->n;

    double A = gsl_vector_get (x, 0);
    double lambda = gsl_vector_get (x, 1);

    size_t i;

    for (i = 0; i < n; i++)
    {
        /* Jacobian matrix J(i,j) = dfi / dxj, */
        /* where fi = (Yi - yi)/sigma[i],      */
        /*       Yi = A * exp(-lambda * i) + b  */
        /* and the xj are the parameters (A,lambda,b) */
        double t = i;
        double e = exp(-lambda * t);
        gsl_matrix_set (J, i, 0, e);
        gsl_matrix_set (J, i, 1, -t * A * e);
        gsl_matrix_set (J, i, 2, 1.0);
    }
    return GSL_SUCCESS;
}

void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
    //return;
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    gsl_vector *x = gsl_multifit_nlinear_position(w);

    double rcond = 1.0;
    /* compute reciprocal condition number of J(x) */
    gsl_multifit_nlinear_rcond(&rcond, w);

    fprintf(stderr, "iter %2zu: l1 = %.4f, l2 = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
            iter,
            gsl_vector_get(x, 0),
            gsl_vector_get(x, 1),
            1.0 / rcond,
            gsl_blas_dnrm2(f));
}

#if 0
TEST(GSL, experiment) {

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();


    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params =
            gsl_multifit_nlinear_default_parameters();
    const size_t n = N;
    const size_t p = 3;

    gsl_vector *f;
    gsl_matrix *J;
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    double y[N], weights[N];
    struct data d = { n, y };
    double x_init[3] = { 1.0, 1.0, 0.0 }; /* starting values */
    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    gsl_vector_view wts = gsl_vector_view_array(weights, n);
    gsl_rng * r;
    double chisq, chisq0;
    int status, info;
    size_t i;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 0.0;

    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);

    /* define the function to be minimized */
    fdf.f = expb_f;
    fdf.df = expb_df;   /* set to NULL for finite-difference Jacobian */
    fdf.fvv = NULL;     /* not using geodesic acceleration */
    fdf.n = n;
    fdf.p = p;
    fdf.params = &d;

    /* this is the data to be fitted */

    for (i = 0; i < n; i++)
    {
        double t = i;
        double yi = 1.0 + 5 * exp (-0.1 * t);
        double si = 0.1 * yi;
        double dy = gsl_ran_gaussian(r, si);

        weights[i] = 1.0 / (si * si);
        y[i] = yi + dy;
        //printf ("data: %zu %g %g\n", i, y[i], si);
    };

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

    /* initialize solver with starting point and weights */
    gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

    /* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

    /* solve the system with a maximum of 20 iterations */
    status = gsl_multifit_nlinear_driver(20, xtol, gtol, ftol,
                                         NULL, NULL, &info, w);

    /* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (J, 0.0, covar);

    /* compute final cost */
    gsl_blas_ddot(f, f, &chisq);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) std::sqrt(gsl_matrix_get(covar,i,i))

    fprintf(stderr, "summary from method '%s/%s'\n",
            gsl_multifit_nlinear_name(w),
            gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "number of iterations: %zu\n",
            gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n",
            (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %f\n", std::sqrt(chisq0));
    fprintf(stderr, "final   |f(x)| = %f\n", std::sqrt(chisq));

    {
        double dof = n - p;
        double c = GSL_MAX_DBL(1, std::sqrt(chisq / dof));

        fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

        fprintf (stderr, "A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
        fprintf (stderr, "lambda = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
        fprintf (stderr, "b      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
    }

    fprintf (stderr, "status = %s\n", gsl_strerror (status));

    gsl_multifit_nlinear_free (w);
    gsl_matrix_free (covar);
    gsl_rng_free (r);

    std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();

    std::cout << "GSL time: " << time << " (" << 1.0/time << "/s)" << std::endl;

}
#endif

struct bezierfit_data {
    size_t n;
    double * x;
    double * y;
    double start_x, start_y, end_x, end_y, mid1_x, mid1_y, mid2_x, mid2_y;
};

double evaluateBezier(double t, double start, double end, double mid1, double mid2, double l1, double l2) {
    return (1-t)*(1-t)*(1-t)*start
            + 3*(1-t)*(1-t)*t*(start+mid1*l1)
            + 3*(1-t)*t*t*(end + mid2*l2)
            + t*t*t*end;
}



int
bezierfit_f (const gsl_vector * x, void *data,
             gsl_vector * f)
{
    struct bezierfit_data * d = static_cast<struct bezierfit_data *>(data);
    size_t n = d->n;

    double l1 = gsl_vector_get (x, 0);
    double l2 = gsl_vector_get (x, 1);

    size_t ii;
    for (ii = 0; ii < n; ii++)
    {
        /* Model Yi = A * exp(-lambda * i) + b */
        /*
        double t = i;
        double Yi = A * exp (-lambda * t) + b;
        gsl_vector_set (f, i, Yi - y[i]);
        */
        double t = gsl_vector_get(x, 2 + ii);
        double c0 = (1-t)*(1-t)*(1-t);
        double c1 = 3*(1-t)*(1-t)*t;
        double c2 = 3*(1-t)*t*t;
        double c3 = t*t*t;

        //*
        double dx_t = - d->x[ii]
                + c0 * d->start_x
                + c1 * (d->start_x + l1 * d->mid1_x)
                + c2 * (d->end_x + l2 * d->mid2_x)
                + c3 * d->end_x;
                //*/
        //double dx_t = -d->x[ii] + evaluateBezier(t, d->start_x, d->end_x, d->mid1_x, d->mid2_x, l1, l2);
        gsl_vector_set(f, 2*ii, dx_t);

        double dy_t = - d->y[ii]
                + c0 * d->start_y
                + c1 * (d->start_y + l1 * d->mid1_y)
                + c2 * (d->end_y + l2 * d->mid2_y)
                + c3 * d->end_y;
        gsl_vector_set(f, 2*ii+1, dy_t);
    }

    return GSL_SUCCESS;
}



TEST(GSL, bezierExperiment) {

    std::chrono::high_resolution_clock::time_point start_t = std::chrono::high_resolution_clock::now();

    Geom::Point start(0,0);
    Geom::Point end(1,0);
    Geom::Point dir1(0,1);
    Geom::Point dir2(0,1);
    Geom::CubicBezier orig(start, start+dir1, end+dir2, end);

    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params =
            gsl_multifit_nlinear_default_parameters();
    const size_t num_points = 25;
    const size_t num_conditions = 2 * num_points;
    const size_t num_params = 2 + num_points;

    gsl_vector *f;
    gsl_matrix *J;
    gsl_matrix *covar = gsl_matrix_alloc (num_params, num_params);
    double target_y[num_points], target_x[num_points], weights[num_conditions];
    double start_x, start_y, end_x, end_y, mid1_x, mid1_y, mid2_x, mid2_y;
    start_x = start.x();
    start_y = start.y();
    end_x = end.x();
    end_y = end.y();
    mid1_x = dir1.x()/4;
    mid1_y = dir1.y()/4;
    mid2_x = dir2.x()/4;
    mid2_y = dir2.y()/4;
    /* this is the data to be fitted */
    size_t ii;
#define BEZ_THRESHOLD (1e-8)
    for (ii = 0; ii < num_points; ii++)
    {
        double t = static_cast<double>(ii) / (num_points - 1);
        Geom::Point target = orig.pointAt(t);

        weights[2 * ii] = 1.0;
        weights[2 * ii + 1] = 1.0;
        target_x[ii] = target.x();
        target_y[ii] = target.y();
        printf ("data: %zu %g %g\n", ii, target_x[ii], target_y[ii]);
    };


    struct bezierfit_data d = { num_points, target_x, target_y, start_x, start_y, end_x, end_y, mid1_x, mid1_y, mid2_x, mid2_y};
    double param_init[num_params];
    param_init[0] = 0;
    param_init[1] = 0;
    for (size_t ii = 0; ii < num_points; ++ii) {
        param_init[ii+2] = static_cast<double>(ii+1) / (num_points + 2);
    }
    gsl_vector_view x = gsl_vector_view_array (param_init, num_params);
    gsl_vector_view wts = gsl_vector_view_array(weights, num_conditions);
    gsl_rng * r;
    double chisq, chisq0;
    int status, info;


    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 0.0;

    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);

    /* define the function to be minimized */
    fdf.f = bezierfit_f;
    fdf.df = NULL;   /* set to NULL for finite-difference Jacobian */
    fdf.fvv = NULL;     /* not using geodesic acceleration */
    fdf.n = num_conditions;
    fdf.p = num_params;
    fdf.params = &d;



    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, num_conditions, num_params);

    /* initialize solver with starting point and weights */
    gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

    /* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

    /* solve the system with a maximum of 20 iterations */
    status = gsl_multifit_nlinear_driver(20, xtol, gtol, ftol,
                                         callback, NULL, &info, w);

    /* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (J, 0.0, covar);

    /* compute final cost */
    gsl_blas_ddot(f, f, &chisq);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) std::sqrt(gsl_matrix_get(covar,i,i))

    fprintf(stderr, "summary from method '%s/%s'\n",
            gsl_multifit_nlinear_name(w),
            gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "number of iterations: %zu\n",
            gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n",
            (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %f\n", std::sqrt(chisq0));
    fprintf(stderr, "final   |f(x)| = %f\n", std::sqrt(chisq));

    {
        double dof = num_conditions - num_params;
        double c = GSL_MAX_DBL(1, std::sqrt(chisq / dof));

        fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

        fprintf (stderr, "l1 = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
        fprintf (stderr, "l2 = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
        for (size_t ii = 0; ii < num_points; ++ii) {
            fprintf (stderr, "t_%u = %.5f +/- %.5f\n", ii, FIT(2+ii), c*ERR(2+ii));
        }
    }

    fprintf (stderr, "status = %s\n", gsl_strerror (status));

    gsl_multifit_nlinear_free (w);
    gsl_matrix_free (covar);
    gsl_rng_free (r);

    std::chrono::high_resolution_clock::time_point stop_t = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double>>(stop_t - start_t).count();

    std::cout << "GSL time: " << time << " (" << 1.0/time << "/s)" << std::endl;

}

struct bezierfit_geom_data {
    size_t n;
    Geom::Point * target;
    Geom::Point start, end, dir1, dir2;
};

int
bezierfit_geom_f (const gsl_vector * x, void *data,
             gsl_vector * f)
{
    struct bezierfit_geom_data * d = static_cast<struct bezierfit_geom_data *>(data);
    size_t n = d->n;

    double l1 = gsl_vector_get (x, 0);
    double l2 = gsl_vector_get (x, 1);
    Geom::CubicBezier curve(d->start, d->start + l1 * d->dir1, d->end + l2 * d->dir2, d->end);

    size_t ii;
    for (ii = 0; ii < n; ii++)
    {
        /* Model Yi = A * exp(-lambda * i) + b */
        /*
        double t = i;
        double Yi = A * exp (-lambda * t) + b;
        gsl_vector_set (f, i, Yi - y[i]);
        */
        double t = gsl_vector_get(x, 2 + ii);

        Geom::Point f_t = curve.pointAt(t);
        Geom::Point diff_t = f_t - d->target[ii];

        gsl_vector_set(f, 2*ii, diff_t.x());

        gsl_vector_set(f, 2*ii+1, diff_t.y());
    }

    return GSL_SUCCESS;
}

TEST(GSL, bezierGeomExperiment) {

    std::chrono::high_resolution_clock::time_point start_t = std::chrono::high_resolution_clock::now();

    Geom::Point start(0,0);
    Geom::Point end(1,0);
    Geom::Point dir1(0,1);
    Geom::Point dir2(0,1);
    Geom::CubicBezier orig(start, start+dir1, end+dir2, end);

    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params =
            gsl_multifit_nlinear_default_parameters();
    const size_t num_points = 25;
    const size_t num_conditions = 2 * num_points;
    const size_t num_params = 2 + num_points;

    gsl_vector *f;
    gsl_matrix *J;
    gsl_matrix *covar = gsl_matrix_alloc (num_params, num_params);
    double weights[num_conditions];
    Geom::Point target[num_points];

    /* this is the data to be fitted */
    size_t ii;
#define BEZ_THRESHOLD (1e-8)
    for (ii = 0; ii < num_points; ii++)
    {
        double t = static_cast<double>(ii) / (num_points - 1);
        weights[2 * ii] = 1.0;
        weights[2 * ii + 1] = 1.0;
        target[ii] = orig.pointAt(t);
        printf ("data: %zu %g %g\n", ii, target[ii].x(), target[ii].y());
    };


    struct bezierfit_geom_data d = { num_points, target, start, end, dir1, dir2};
    double param_init[num_params];
    param_init[0] = 0;
    param_init[1] = 0;
    for (size_t ii = 0; ii < num_points; ++ii) {
        param_init[ii+2] = static_cast<double>(ii+1) / (num_points + 2);
    }
    gsl_vector_view x = gsl_vector_view_array (param_init, num_params);
    gsl_vector_view wts = gsl_vector_view_array(weights, num_conditions);
    gsl_rng * r;
    double chisq, chisq0;
    int status, info;


    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 0.0;

    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);

    /* define the function to be minimized */
    fdf.f = bezierfit_geom_f;
    fdf.df = NULL;   /* set to NULL for finite-difference Jacobian */
    fdf.fvv = NULL;     /* not using geodesic acceleration */
    fdf.n = num_conditions;
    fdf.p = num_params;
    fdf.params = &d;



    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, num_conditions, num_params);

    /* initialize solver with starting point and weights */
    gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

    /* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

    /* solve the system with a maximum of 20 iterations */
    status = gsl_multifit_nlinear_driver(20, xtol, gtol, ftol,
                                         callback, NULL, &info, w);

    /* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (J, 0.0, covar);

    /* compute final cost */
    gsl_blas_ddot(f, f, &chisq);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) std::sqrt(gsl_matrix_get(covar,i,i))

    fprintf(stderr, "summary from method '%s/%s'\n",
            gsl_multifit_nlinear_name(w),
            gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "number of iterations: %zu\n",
            gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n",
            (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %f\n", std::sqrt(chisq0));
    fprintf(stderr, "final   |f(x)| = %f\n", std::sqrt(chisq));

    {
        double dof = num_conditions - num_params;
        double c = GSL_MAX_DBL(1, std::sqrt(chisq / dof));

        fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

        fprintf (stderr, "l1 = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
        fprintf (stderr, "l2 = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
        for (size_t ii = 0; ii < num_points; ++ii) {
            fprintf (stderr, "t_%u = %.5f +/- %.5f\n", ii, FIT(2+ii), c*ERR(2+ii));
        }
    }

    fprintf (stderr, "status = %s\n", gsl_strerror (status));

    gsl_multifit_nlinear_free (w);
    gsl_matrix_free (covar);
    gsl_rng_free (r);

    std::chrono::high_resolution_clock::time_point stop_t = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double>>(stop_t - start_t).count();

    std::cout << "GSL time: " << time << " (" << 1.0/time << "/s)" << std::endl;

}

struct BezierfitGeomDistanceData {
    size_t n;
    double width;
    Geom::Point * target;
    Geom::Point start, end, dir1, dir2;
};

int
BezierfitGeomDistanceF (const gsl_vector * x, void *data,
             gsl_vector * f)
{
    struct BezierfitGeomDistanceData * d = static_cast<struct BezierfitGeomDistanceData *>(data);
    size_t const n = d->n;

    double const l1 = gsl_vector_get (x, 0);
    double const l2 = gsl_vector_get (x, 1);
    Geom::CubicBezier const curve(d->start, d->start + l1 * d->dir1, d->end + l2 * d->dir2, d->end);

    size_t ii;
    for (ii = 0; ii < n; ii++)
    {
        /* Model Yi = A * exp(-lambda * i) + b */
        /*
        double t = i;
        double Yi = A * exp (-lambda * t) + b;
        gsl_vector_set (f, i, Yi - y[i]);
        */
        double t = gsl_vector_get(x, 2 + ii);
        if (t < 0) {
            t = 0;
        }
        else if (t > 1) {
            t = 1;
        }

        std::vector<Geom::Point> f_t = curve.pointAndDerivatives(t, 1);

        Geom::Point diff_t = f_t[0] - d->target[ii];

        gsl_vector_set(f, 2*ii, diff_t.length() - d->width);

        double const angle = Geom::dot(f_t[1], diff_t) / (f_t[1].length() * diff_t.length());
        gsl_vector_set(f, 2*ii+1, angle);
    }

    return GSL_SUCCESS;
}

#define GSL_DBG_MSG 0

TEST(GSL, bezierGeomDistanceExperiment) {

    std::chrono::high_resolution_clock::time_point start_t = std::chrono::high_resolution_clock::now();

    Geom::Point start(0,0);
    Geom::Point end(1,0);
    Geom::Point dir1(0,1);
    Geom::Point dir2(0,1);

    Geom::Point offset_start(-1,0);
    Geom::Point offset_end(2,0);
    double const width = 1.0;
    Geom::CubicBezier orig(start, start+dir1, end+dir2, end);

    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params =
            gsl_multifit_nlinear_default_parameters();
    const size_t num_points = 3;
    const size_t num_conditions = 2 * num_points;
    const size_t num_params = 2 + num_points;

    gsl_vector *f;
    gsl_matrix *J;
    gsl_matrix *covar = gsl_matrix_alloc (num_params, num_params);
    double weights[num_conditions];
    Geom::Point target[num_points];

    /* this is the data to be fitted */
    size_t ii;
#define BEZ_THRESHOLD (1e-8)
    for (ii = 0; ii < num_points; ii++)
    {
        double t = static_cast<double>(ii+1) / (num_points +2);
        weights[2 * ii] = 1;
        weights[2 * ii + 1] = 1;
        target[ii] = orig.pointAt(t);
#if GSL_DBG_MSG
        printf ("data: %zu %g %g\n", ii, target[ii].x(), target[ii].y());
#endif
    };



    struct BezierfitGeomDistanceData d = { num_points, width, target, offset_start, offset_end, dir1, dir2};
    double param_init[num_params];
    param_init[0] = 1;
    param_init[1] = 1;
    for (size_t ii = 0; ii < num_points; ++ii) {
        param_init[ii+2] = static_cast<double>(ii+1) / (num_points +2);
    }
    gsl_vector_view x = gsl_vector_view_array (param_init, num_params);
    gsl_vector_view wts = gsl_vector_view_array(weights, num_conditions);
    gsl_rng * r;
    double chisq, chisq0;
    int status, info;


    const double xtol = 1e-4;
    const double gtol = 1e-4;
    const double ftol = 0.0;

    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);

    /* define the function to be minimized */
    fdf.f = BezierfitGeomDistanceF;
    fdf.df = NULL;   /* set to NULL for finite-difference Jacobian */
    fdf.fvv = NULL;     /* not using geodesic acceleration */
    fdf.n = num_conditions;
    fdf.p = num_params;
    fdf.params = &d;



    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, num_conditions, num_params);

    /* initialize solver with starting point and weights */
    gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

#if GSL_DBG_MSG
    /* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);
#endif

    /* solve the system with a maximum of 30 iterations */
    status = gsl_multifit_nlinear_driver(30, xtol, gtol, ftol,
                                         NULL, NULL, &info, w);


#if GSL_DBG_MSG
    /* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (J, 0.0, covar);

    /* compute final cost */
    gsl_blas_ddot(f, f, &chisq);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) std::sqrt(gsl_matrix_get(covar,i,i))
    fprintf(stderr, "summary from method '%s/%s'\n",
            gsl_multifit_nlinear_name(w),
            gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "number of iterations: %zu\n",
            gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n",
            (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %f\n", std::sqrt(chisq0));
    fprintf(stderr, "final   |f(x)| = %f\n", std::sqrt(chisq));

    {
        double dof = num_conditions - num_params;
        double c = GSL_MAX_DBL(1, std::sqrt(chisq / dof));

        fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

        fprintf (stderr, "l1 = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
        fprintf (stderr, "l2 = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
        for (size_t ii = 0; ii < num_points; ++ii) {
            fprintf (stderr, "t_%u = %.5f +/- %.5f\n", ii, FIT(2+ii), c*ERR(2+ii));
        }
    }

    fprintf (stderr, "status = %s\n", gsl_strerror (status));
#endif
    gsl_multifit_nlinear_free (w);
    gsl_matrix_free (covar);
    gsl_rng_free (r);

    std::chrono::high_resolution_clock::time_point stop_t = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double>>(stop_t - start_t).count();

    std::cout << "GSL time: " << time << " (" << 1.0/time << "/s)" << std::endl;

}

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . +))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=99 :
