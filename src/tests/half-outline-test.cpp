/** @file
 * @brief Unit tests for Affine
 * Uses the Google Testing Framework
 *//*
 * Authors:
 *   Krzysztof Kosiński <tweenk.pl@gmail.com>
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

#include <gtest/gtest.h>
#include <2geom/affine.h>
#include <2geom/transforms.h>
#include <2geom/svg-path-parser.h>
#include <2geom/svg-path-writer.h>
#include <2geom/sbasis-to-bezier.h>

#include "libRunningStats/runningstats.h"
#include <boost/filesystem.hpp>
#include "../toys/helper/geom-pathstroke.h"

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
        {
            Path new_sub;
            Inkscape::offset_curve(new_sub, &orig[ii], width);
            Path orig_sub;
            orig_sub.append(orig[ii]);
            auto stats = symmetricDistanceStats(new_sub, orig_sub, width, num_tests);
            std::cout << "subpath: " << ii << ", forward: " << stats.print() << std::endl;
            std::cout << "New data: " << write_svg_path(new_sub) << std::endl;
        }
        {
            Path new_sub;
            Inkscape::offset_curve(new_sub, orig[ii].reverse(), width);
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
