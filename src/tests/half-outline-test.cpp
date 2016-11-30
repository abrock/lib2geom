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

#define NUM_SPEED_TESTS 50

void speedTest(std::string const dirName) {
    double const width = 1;
    double const miter = 0;
    const auto files = list_files(dirName);
    std::vector<Geom::Path> paths;
    for (const auto& file : files) {
        if (file.extension() == ".svgd") {
            Geom::PathVector tmp = read_svgd(file.c_str());
            ASSERT_EQ(1, tmp.size());
            paths.push_back(tmp[0]);
        }
    }
    std::cout << "Running speedTest with file in " << dirName << std::endl;
    std::cout << "New method: ";
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    size_t counter = 0;
    for (size_t ii = 0; ii < NUM_SPEED_TESTS; ++ii) {
        for (auto& path: paths) {
            auto result = Inkscape::half_outline(path, width, miter);
            counter += path.size();
        }
    }
    std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    double const new_speed = static_cast<double>(NUM_SPEED_TESTS * paths.size()) / time;
    double const new_speed_curves = static_cast<double>(NUM_SPEED_TESTS * counter) / time;
    std::cout << new_speed << " paths per second (" << new_speed_curves << " curves/s)" << std::endl;


    std::cout << "Old method: ";
    start = std::chrono::high_resolution_clock::now();
    counter = 0;
    for (size_t ii = 0; ii < NUM_SPEED_TESTS; ++ii) {
        for (auto& path: paths) {
            auto result = half_outline_old(path, width, miter);
            counter += path.size();
        }
    }
    stop = std::chrono::high_resolution_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    double const old_speed = static_cast<double>(NUM_SPEED_TESTS * paths.size()) / time;
    double const old_speed_curves = static_cast<double>(NUM_SPEED_TESTS * counter) / time;
    std::cout << old_speed << " paths per second (" << old_speed_curves << " curves/s)" << std::endl;

    std::cout << "Old method is " << old_speed / new_speed << " times as fast as new method" << std::endl;
    std::cout << std::endl;
}

TEST(HalfOutlineTest, speedTest) {
    speedTest("half-outline-degenerate-paths");
    speedTest("half-outline-difficult-paths");
    speedTest("half-outline-embroidery-paths");
    speedTest("half-outline-known-issues");
    speedTest("half-outline-simple-paths");
    speedTest("half-outline-many-nodes");
}

void speedTestSubCurves(std::string const dirName) {
    return;

    double const width = 1;
    const auto files = list_files(dirName);
    std::vector<Geom::Path> paths;
    for (const auto& file : files) {
        if (file.extension() == ".svgd") {
            Geom::PathVector tmp = read_svgd(file.c_str());
            ASSERT_EQ(1, tmp.size());
            paths.push_back(tmp[0]);
        }
    }
    std::cout << "Running speedTest with file in " << dirName << std::endl;
    std::cout << "New method: ";
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    size_t counter = 0;
    for (size_t ii = 0; ii < NUM_SPEED_TESTS; ++ii) {
        for (auto& path: paths) {
            for (size_t ii = 0; ii < path.size(); ++ii) {
                Geom::Path tmp;
                Inkscape::offset_curve(tmp, &path[ii], width);
                counter++;
            }
        }
    }
    std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    double const new_speed = static_cast<double>(counter) / time;
    std::cout << new_speed << " curves per second" << std::endl;


    std::cout << "Old method: ";
    start = std::chrono::high_resolution_clock::now();
    counter = 0;
    for (size_t ii = 0; ii < NUM_SPEED_TESTS; ++ii) {
        for (auto& path: paths) {
            for (size_t ii = 0; ii < path.size(); ++ii) {
                Geom::Path tmp;
                offset_curve_old(tmp, &path[ii], width);
                counter++;
            }
        }
    }
    stop = std::chrono::high_resolution_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    double const old_speed = static_cast<double>(counter) / time;
    std::cout << old_speed << " curves per second" << std::endl;

    std::cout << "Old method is " << old_speed / new_speed << " times as fast as new method" << std::endl;
    std::cout << std::endl;
}

TEST(HalfOutlineTest, speedTestSubCurves) {
    speedTestSubCurves("half-outline-degenerate-paths");
    speedTestSubCurves("half-outline-difficult-paths");
    speedTestSubCurves("half-outline-embroidery-paths");
    speedTestSubCurves("half-outline-known-issues");
    speedTestSubCurves("half-outline-simple-paths");
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
