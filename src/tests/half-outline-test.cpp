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

namespace Geom {

#include "half-outline-old.cpp"

namespace fs = boost::filesystem;

std::vector<fs::path> list_files(const fs::path & dir_path) {
    std::vector<fs::path> result;
    if (!exists(dir_path)) {
        return result;
    }
    fs::directory_iterator end_itr; // default construction yields past-the-end
    for (fs::directory_iterator itr(dir_path);
         itr != end_itr;
         ++itr) {
        if (fs::is_regular_file(itr->status())) {
            result.push_back(itr->path());
        }
    }
    return result;
}

/**
 * @brief distanceStats calculates statistics of distances between
 * @param result
 * @param a
 * @param b
 * @param expected_width
 * @param num_test_points
 */
void distanceStats(
        QuantileStats<double>& result,
        Geom::Path const & a,
        Geom::Path const & b,
        double expected_width,
        size_t const num_test_points = 1024)
{
    expected_width = std::abs(expected_width);
    for (size_t jj = 0; jj < a.size(); ++jj) {
        for (size_t ii = 0; ii <= num_test_points; ++ii) {
            Geom::PathTime const t(jj, static_cast<double>(ii) / num_test_points);
            Geom::Point const a_point = a.pointAt(t);
            Geom::Point const b_point = b.pointAt(b.nearestTime(a_point));
            double const diff = (a_point - b_point).length();
            result.push(100 * (diff - expected_width)/expected_width);
        }
    }
}

QuantileStats<double> symmetricDistanceStats(
        Geom::Path const & a,
        Geom::Path const & b,
        double const expected_width,
        size_t const num_test_points = 1024) {
    QuantileStats<double> result;
    distanceStats(result, a, b, expected_width, num_test_points);
    distanceStats(result, b, a, expected_width, num_test_points);
    return result;
}

void halfOutlineTestNoExplosion(const fs::path& name, const double width = 5, const double miter = 0) {
    Geom::PathVector const orig_pathvector = read_svgd(name.c_str());
    ASSERT_EQ(1u, orig_pathvector.size());
    Geom::Path const orig = orig_pathvector[0];
    Geom::Path const result = Inkscape::half_outline(orig, width, miter);
    auto stats = symmetricDistanceStats(orig, result, width);
    EXPECT_NEAR(stats.getMin(), 0, width);
    EXPECT_NEAR(stats.getMax(), 0, width);

#define DEBUG_VERBOSE 1
#if (DEBUG_VERBOSE)
    std::cout << "size: " << orig.size() << " name: " << name << std::endl;
    std::cout << "result size: " << result.size() << std::endl;
    std::cout << "stats: " << stats.print()
              << ", 20/80: " << stats.getQuantile(.2) << ", " << stats.getQuantile(.8) << std::endl;
#endif


    Geom::Path const result_old = half_outline_old(orig, width, miter);
    stats = symmetricDistanceStats(orig, result_old, width);
    EXPECT_NEAR(stats.getMin(), 0, width);
    EXPECT_NEAR(stats.getMax(), 0, width);

#if (DEBUG_VERBOSE)
    std::cout << "old result size: " << result_old.size() << std::endl;
    std::cout << "stats: " << stats.print()
              << ", 20/80: " << stats.getQuantile(.2) << ", " << stats.getQuantile(.8) << std::endl;
#endif

    Geom::Path const result_rev = Inkscape::half_outline(orig.reversed(), width, miter);
    stats = symmetricDistanceStats(orig, result_rev, width);
    EXPECT_NEAR(stats.getMin(), 0, width);
    EXPECT_NEAR(stats.getMax(), 0, width);

#if (DEBUG_VERBOSE)
    std::cout << "rev. size: " << result_rev.size() << std::endl;
    std::cout << "rev. stats: " << stats.print()
              << ", 20/80: " << stats.getQuantile(.2) << ", " << stats.getQuantile(.8) << std::endl;
#endif

    Geom::Path const result_rev_old = half_outline_old(orig.reversed(), width, miter);
    stats = symmetricDistanceStats(orig, result_rev_old, width);
    EXPECT_NEAR(stats.getMin(), 0, width);
    EXPECT_NEAR(stats.getMax(), 0, width);

#if (DEBUG_VERBOSE)
    std::cout << "old rev. size: " << result_rev_old.size() << std::endl;
    std::cout << "old rev. stats: " << stats.print()
              << ", 20/80: " << stats.getQuantile(.2) << ", " << stats.getQuantile(.8) << std::endl;

    std::cout << std::endl;
#endif
}

void offsetCurveTest(const fs::path& name, const double width = 5, const double miter = 0) {
    Geom::PathVector const orig_pathvector = read_svgd(name.c_str());
    ASSERT_EQ(1u, orig_pathvector.size());
    Geom::Path const orig = orig_pathvector[0];

    const size_t num_tests = 1024;

    std::cout << "file: " << name << std::endl;
    for (size_t ii = 0; ii < orig.size(); ++ii) {
        Path orig_path;
        orig_path.append(orig[ii]);
        std::cout << "Data: " << write_svg_path(orig_path) << std::endl;
        double const tolerance = 5.0 * (width / 100);
        {
            Path new_sub;
            Inkscape::offset_curve(new_sub, &orig[ii], width, tolerance);
            Path orig_sub;
            orig_sub.append(orig[ii]);
            auto stats = symmetricDistanceStats(new_sub, orig_sub, width, num_tests);
            std::cout << "subpath: " << ii << ", forward: " << stats.print() << std::endl;
            std::cout << "New data: " << write_svg_path(new_sub) << std::endl;
        }
        {
            Path new_sub;
            Inkscape::offset_curve(new_sub, orig[ii].reverse(), width, tolerance);
            Path orig_sub;
            orig_sub.append(orig[ii]);
            auto stats = symmetricDistanceStats(new_sub, orig_sub, width, num_tests);
            std::cout << "subpath: " << ii << ", reverse: " << stats.print() << std::endl;
            std::cout << "New data: " << write_svg_path(new_sub) << std::endl;
        }
        std::cout << std::endl;
    }

}


TEST(HalfOutlineTest, noExplosion) {
    const auto files = list_files("half-outline-difficult-paths");
    for (const auto& file : files) {
        halfOutlineTestNoExplosion(file);
    }
}

TEST(HalfOutlineTest, knownIssues) {
    const auto files = list_files("half-outline-known-issues");
    for (const auto& file : files) {
        offsetCurveTest(file, .1);
    }
}

#define COMBINATORICAL_EXPLOSION_ERROR_LIMIT 20
#define COMBINATORICAL_EXPLOSIONS_MAXIT 400
#define COMBINATORICAL_EXPLOSIONS_MAX_NODEFACTOR 40

int findCombinatoricalExplosion(
        Geom::Path const input,
        std::ostream& out,
        double const width = .5,
        double const miter = 0,
        bool const use_old_algorithm = false) {
    Geom::Path last_path = input;
    out << "<g>";
    for (size_t ii = 0; ii < COMBINATORICAL_EXPLOSIONS_MAXIT; ++ii) {
        Geom::Path result;
        if (use_old_algorithm) {
            result = half_outline_old(last_path, width, miter);
        }
        else {
            result = Inkscape::half_outline(last_path, width, miter);
        }
        out << "<path style=\"display:inline;fill:none;fill-opacity:1;stroke:#aaaaaa;stroke-width:" << width/2 << ";stroke-miterlimit:6;stroke-dasharray:none;stroke-opacity:1;enable-background:new\" d=\""
            << write_svg_path(result) << "\"/>" << std::endl;
        auto stats = symmetricDistanceStats(last_path, result, width, 50);
        //EXPECT_NEAR(stats.getMin(), 0, 1);
        //EXPECT_NEAR(stats.getMax(), 0, 1);
        //std::cout << ii << "\t" << result.size() << std::endl;

        if (result.size() > COMBINATORICAL_EXPLOSIONS_MAX_NODEFACTOR * input.size()) {
            return ii;
        }
        if (std::abs(stats.getMin()) > COMBINATORICAL_EXPLOSION_ERROR_LIMIT) {
            return ii;
        }
        if (std::abs(stats.getMax()) > COMBINATORICAL_EXPLOSION_ERROR_LIMIT) {
            return ii;
        }
        last_path = result;
    }
    out << "</g>";
    return -1;
}

void findCombinatoricalExplosion(
        Geom::Path const input,
        std::string filename,
        double const width = .5,
        double const miter = 0,
        bool const use_old_algorithm = false) {

    std::ofstream out(filename);
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><svg "
        << "xmlns:dc=\"http://purl.org/dc/elements/1.1/\"  xmlns:cc=\"http://creativecommons.org/ns#\"  xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"  xmlns:svg=\"http://www.w3.org/2000/svg\"  xmlns=\"http://www.w3.org/2000/svg\"  xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\"  xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\"  width=\"400mm\"  height=\"350mm\"  viewBox=\"0 0 400 350\"  id=\"svg12765\"  version=\"1.1\">"
        << std::endl;
    std::cout << "Forward explosion: " << findCombinatoricalExplosion(input, out, width, miter, use_old_algorithm) << std::endl;
    std::cout << "Reverse explosion: " << findCombinatoricalExplosion(input.reversed(), out, width, miter, use_old_algorithm) << std::endl;

    out << "</svg>";
}

void repeatedOffsetTest(const fs::path& name, const double width = 5, const double miter = 0) {
    Geom::PathVector const orig_pathvector = read_svgd(name.c_str());
    ASSERT_EQ(1u, orig_pathvector.size());
    Geom::Path const orig = orig_pathvector[0];

    std::cout << "File: " << name << std::endl;
    std::cout << "New algorithm:" << std::endl;
    findCombinatoricalExplosion(orig, name.string() + "-new.svg", width, miter, false);
    std::cout << "Old algorithm:" << std::endl;
    findCombinatoricalExplosion(orig, name.string() + "-old.svg", width, miter, true);
    std::cout << std::endl;
}

TEST(HalfOutlineTest, combinatoricalExplosion) {
    const auto files = list_files("half-outline-embroidery-paths");
    for (const auto& file : files) {
        if (file.extension() == ".svgd") {
            repeatedOffsetTest(file, .5);
        }
    }
}


#define NUM_SPEED_TESTS 50

TEST(HalfOutlineTest, speedTest) {
    double const width = 1;
    double const miter = 0;
    const auto files = list_files("half-outline-embroidery-paths");
    std::vector<Geom::Path> paths;
    for (const auto& file : files) {
        if (file.extension() == ".svgd") {
            Geom::PathVector tmp = read_svgd(file.c_str());
            ASSERT_EQ(1, tmp.size());
            paths.push_back(tmp[0]);
        }
    }
    std::cout << "New method:" << std::endl;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (size_t ii = 0; ii < NUM_SPEED_TESTS; ++ii) {
        for (auto& path: paths) {
            auto result = Inkscape::half_outline(path, width, miter);
        }
    }
    std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    double const new_speed = static_cast<double>(NUM_SPEED_TESTS * paths.size()) / time;
    std::cout << new_speed << " paths per second" << std::endl;


    std::cout << "Old method:" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    for (size_t ii = 0; ii < NUM_SPEED_TESTS; ++ii) {
        for (auto& path: paths) {
            auto result = half_outline_old(path, width, miter);
        }
    }
    stop = std::chrono::high_resolution_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    double const old_speed = static_cast<double>(NUM_SPEED_TESTS * paths.size()) / time;
    std::cout << old_speed << " paths per second" << std::endl;

    std::cout << "Old method is " << old_speed / new_speed << " times as fast as new method" << std::endl;
}


TEST(HalfOutlineTest, precisionComparison) {
    double const width = 1;
    double const miter = 0;
    const auto files = list_files("half-outline-embroidery-paths");
    for (const auto& file : files) {
        if (file.extension() == ".svgd") {
            Geom::PathVector tmp = read_svgd(file.c_str());
            ASSERT_EQ(1, tmp.size());
            auto path = tmp[0];
            std::cout << "File: " << file << std::endl;
            {
                auto result = Inkscape::half_outline(path, width, miter);
                auto stats = symmetricDistanceStats(result, path, width);
                std::cout << "New: " << result.size() << "\t" << stats.print() << std::endl;
            }
            {
                auto result = half_outline_old(path, width, miter);
                auto stats = symmetricDistanceStats(result, path, width);
                std::cout << "Old: " << result.size() << "\t" << stats.print() << std::endl;
            }
            std::cout << std::endl;
        }
    }
}

TEST(HalfOutlineTest, degeneratePaths) {
    double const width = .1;
    double const miter = 0;
    const auto files = list_files("half-outline-degenerate-paths");
    for (const auto& file : files) {
        if (file.extension() == ".svgd") {
            Geom::PathVector tmp = read_svgd(file.c_str());
            ASSERT_EQ(1, tmp.size());
            auto const path = tmp[0];
            std::cout << "File: " << file << ": " << path.size() << std::endl;
            //std::cout << write_svg_path(path) << std::endl;
            std::cout << "Forward tests:" << std::endl;
            {
                Geom::Path result = Inkscape::half_outline(path, width, miter);
                //std::cout << write_svg_path(result) << std::endl;
                auto stats = symmetricDistanceStats(result, path, width);
                std::cout << "New: " << result.size() << "\t" << stats.print() << std::endl;
            }
            {
                Geom::Path result = half_outline_old(path, width, miter);
                //std::cout << write_svg_path(result) << std::endl;
                auto stats = symmetricDistanceStats(result, path, width);
                std::cout << "Old: " << result.size() << "\t" << stats.print() << std::endl;
            }
            std::cout << "Reversed tests:" << std::endl;
            {
                Geom::Path result = Inkscape::half_outline(path.reversed(), width, miter);
                //std::cout << write_svg_path(result) << std::endl;
                auto stats = symmetricDistanceStats(result, path, width);
                std::cout << "New: " << result.size() << "\t" << stats.print() << std::endl;
            }
            {
                Geom::Path result = half_outline_old(path.reversed(), width, miter);
                //std::cout << write_svg_path(result) << std::endl;
                auto stats = symmetricDistanceStats(result, path, width);
                std::cout << "Old: " << result.size() << "\t" << stats.print() << std::endl;
            }
            std::cout << std::endl;
        }
    }
}



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
    return;
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    gsl_vector *x = gsl_multifit_nlinear_position(w);
    double rcond = 1.0;

    /* compute reciprocal condition number of J(x) */
    //gsl_multifit_nlinear_rcond(&rcond, w);

    fprintf(stderr, "iter %2zu: A = %.4f, lambda = %.4f, b = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
            iter,
            gsl_vector_get(x, 0),
            gsl_vector_get(x, 1),
            gsl_vector_get(x, 2),
            1.0 / rcond,
            gsl_blas_dnrm2(f));
}
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

} // end namespace Geom

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
