/** @file
 * @brief Unit tests for Line and related functions
 * Uses the Google Testing Framework
 *//*
 * Authors:
 *   Alexander Brock <alexander@lunar-orbit.de>
 *
 * Copyright 2015 Authors
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

#include "testing.h"
#include <iostream>
#include <glib.h>

#include <2geom/line.h>
#include <2geom/affine.h>

using namespace Geom;

TEST(NearestTimesTest, Speed) {
    std::vector<std::pair<Path, Path> > testset;
    {
        Path a, b;
        a.append(CubicBezier(Point(0,0), Point(1,0), Point(1,1), Point(0,1)));
        b.append(CubicBezier(Point(1,0), Point(0,0), Point(0,1), Point(1,1)));
        testset.push_back(std::make_pair(a, b));
    }

    for (std::pair<Path, Path> const& test : testset) {
        for (size_t ii = 0; ii < 10000; ++ii) {
            std::pair<PathTime, PathTime> result = nearest_times(test.first, test.second);
        }
    }
}

TEST(NearestTimesTest, TouchingCubicBezierSimple) {
    Path a, b;
    a.append(CubicBezier(Point(0,0), Point(0, 8./6.), Point(1, 8./6.), Point(1,0)));
    b.append(CubicBezier(Point(0,2), Point(0, 4./6.), Point(1, 4./6.), Point(1,2)));
    {
        double dist = std::numeric_limits<Coord>::max();
        std::pair<PathTime, PathTime> result = nearest_times(a,b, &dist);
        EXPECT_NEAR(dist, 0.0, 1e-10);
        EXPECT_NEAR(result.first.t, .5, 1e-10);
        EXPECT_NEAR(result.second.t, .5, 1e-10);
        EXPECT_EQ(result.first.curve_index, 0);
        EXPECT_EQ(result.second.curve_index, 0);
    }
}

TEST(NearestTimesTest, TouchingCubicBezierDivided) {
    Path a, b;
    CubicBezier bez_a(Point(0,0), Point(0, 8./6.), Point(1, 8./6.), Point(1,0));
    CubicBezier bez_b(Point(0,2), Point(0, 4./6.), Point(1, 4./6.), Point(1,2));
    std::pair<CubicBezier, CubicBezier> div_a, div_b;
    div_a = bez_a.subdivide(.3);
    div_b = bez_b.subdivide(.7);

    a.append(div_a.first);
    a.append(div_a.second);

    b.append(div_b.first);
    b.append(div_b.second);

    {
        double dist = std::numeric_limits<Coord>::max();
        std::pair<PathTime, PathTime> result = nearest_times(a,b, &dist);
        EXPECT_NEAR(dist, 0.0, 1e-10);
        EXPECT_EQ(result.first.curve_index, 1);
        EXPECT_EQ(result.second.curve_index, 0);
    }
}

TEST(NearestTimesTest, TouchingQuadraticBezierSimple) {
    {
        Path a, b;
        a.append(QuadraticBezier(Point(0,0), Point(0, 2), Point(1,0)));
        b.append(QuadraticBezier(Point(0,2), Point(0, 0), Point(1,2)));
        double dist = std::numeric_limits<Coord>::max();
        std::pair<PathTime, PathTime> result = nearest_times(a,b, &dist);
        EXPECT_NEAR(dist, 0.0, 1e-10);
        EXPECT_NEAR(result.first.t, .5, 1e-10);
        EXPECT_NEAR(result.second.t, .5, 1e-10);
        EXPECT_EQ(result.first.curve_index, 0);
        EXPECT_EQ(result.second.curve_index, 0);
    }
    {
        Path a, b;
        a.append(QuadraticBezier(Point(0,0), Point(0, 2), Point(0,0)));
        b.append(QuadraticBezier(Point(0,2), Point(0, 0), Point(0,2)));
        double dist = std::numeric_limits<Coord>::max();
        std::pair<PathTime, PathTime> result = nearest_times(a,b, &dist);
        EXPECT_NEAR(dist, 0.0, 1e-10);
        EXPECT_NEAR(result.first.t, .5, 1e-10);
        EXPECT_NEAR(result.second.t, .5, 1e-10);
        EXPECT_EQ(result.first.curve_index, 0);
        EXPECT_EQ(result.second.curve_index, 0);
    }
}


TEST(NearestTimesTest, TouchingQuadraticBezierDivided) {
    Path a, b;
    QuadraticBezier bez_a(Point(0,0), Point(0, 2), Point(1,0));
    QuadraticBezier bez_b(Point(0,2), Point(0, 0), Point(1,2));
    std::pair<QuadraticBezier, QuadraticBezier> div_a, div_b;
    div_a = bez_a.subdivide(.3);
    div_b = bez_b.subdivide(.7);

    a.append(div_a.first);
    a.append(div_a.second);

    b.append(div_b.first);
    b.append(div_b.second);

    {
        double dist = std::numeric_limits<Coord>::max();
        std::pair<PathTime, PathTime> result = nearest_times(a,b, &dist);
        EXPECT_NEAR(dist, 0.0, 1e-4);
        EXPECT_EQ(result.first.curve_index, 1);
        EXPECT_EQ(result.second.curve_index, 0);
    }
}
