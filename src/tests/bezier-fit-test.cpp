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

#include <2geom/bezier-utils.h>

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


TEST(CubicBezier, fitTest) {
    Geom::Point start(0,0);
    Geom::Point mid1(0,1);
    Geom::Point mid2(1,-1);
    Geom::Point end(1,0);
    Geom::CubicBezier bez(start, mid1, mid2, end);

    std::vector<Geom::Point> target;
    const size_t num_points = 20;
    for (size_t ii = 0; ii < num_points; ++ii) {
        double const t = static_cast<double>(ii) / (num_points - 1);
        target.push_back(bez.pointAt(t));
    }
    Geom::CubicBezier fitted(start, start, end, end);
    auto result = Geom::fit_bezier(fitted, target);

    std::cout << "Running speedTest" << std::endl;
    std::cout << "New method: ";
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    size_t counter = 0;
    for (size_t ii = 0; ii < 300; ++ii) {
        Geom::CubicBezier fitted(start, start, end, end);
        result = Geom::fit_bezier(fitted, target);
        counter++;
    }
    std::chrono::high_resolution_clock::time_point stop_time = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double>>(stop_time - start_time).count();
    double const new_speed = static_cast<double>(counter) / time;
    std::cout << new_speed << " curves per second" << std::endl;
    std::cout << "Worst error: " << result.first << " at t=" << result.second << std::endl;
}


TEST(CubicBezier, degenerateCurveTest) {
    Geom::Point start(0,0);
    Geom::Point mid1(-5,0);
    Geom::Point mid2(3, 0);
    Geom::Point end(5,0);
    Geom::CubicBezier bez(start, mid1, mid2, end);

    std::vector<Geom::Point> target;
    const size_t num_points = 20;
    for (size_t ii = 0; ii < num_points; ++ii) {
        double const t = static_cast<double>(ii) / (num_points - 1);
        target.push_back(bez.pointAt(t));
    }
    Geom::CubicBezier fitted(start, start, end, end);
    auto result = Geom::fit_bezier(fitted, target);

    std::cout << "Running speedTest" << std::endl;
    std::cout << "New method: ";
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    size_t counter = 0;
    for (size_t ii = 0; ii < 300; ++ii) {
        Geom::CubicBezier fitted(start, start, end, end);
        result = Geom::fit_bezier(fitted, target);
        counter++;
    }
    std::chrono::high_resolution_clock::time_point stop_time = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double>>(stop_time - start_time).count();
    double const new_speed = static_cast<double>(counter) / time;
    std::cout << new_speed << " curves per second" << std::endl;
    std::cout << "Worst error: " << result.first << " at t=" << result.second << std::endl;
}

TEST(CubicBezier, nearlyDegenerateCurveTest) {
    Geom::Point start(0,0);
    Geom::Point mid1(-5,0.1);
    Geom::Point mid2(3, 0);
    Geom::Point end(5,0);
    Geom::CubicBezier bez(start, mid1, mid2, end);

    std::vector<Geom::Point> target;
    const size_t num_points = 20;
    for (size_t ii = 0; ii < num_points; ++ii) {
        double const t = static_cast<double>(ii) / (num_points - 1);
        target.push_back(bez.pointAt(t));
    }
    Geom::CubicBezier fitted(start, start, end, end);
    auto result = Geom::fit_bezier(fitted, target);

    std::cout << "Running speedTest" << std::endl;
    std::cout << "New method: ";
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    size_t counter = 0;
    for (size_t ii = 0; ii < 300; ++ii) {
        Geom::CubicBezier fitted(start, start, end, end);
        result = Geom::fit_bezier(fitted, target);
        counter++;
    }
    std::chrono::high_resolution_clock::time_point stop_time = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double>>(stop_time - start_time).count();
    double const new_speed = static_cast<double>(counter) / time;
    std::cout << new_speed << " curves per second" << std::endl;
    std::cout << "Worst error: " << result.first << " at t=" << result.second << std::endl;
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
