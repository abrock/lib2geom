/** @file
 * @brief Nearest time routines for D2<SBasis> and Piecewise<D2<SBasis>>
 *//*
 * Authors:
 *   Marco Cecchetti <mrcekets at gmail.com>
 *
 * Copyright 2007-2008  authors
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


#include <2geom/nearest-time.h>
#include <algorithm>
#include <2geom/curve.h>
#include <2geom/path.h>


namespace Geom
{

Coord nearest_time(Point const &p, D2<Bezier> const &input, Coord from, Coord to)
{
    Interval domain(from, to);
    bool partial = false;

    if (domain.min() < 0 || domain.max() > 1) {
        THROW_RANGEERROR("[from,to] interval out of bounds");
    }

    if (input.isConstant(0)) return from;

    D2<Bezier> bez;
    if (domain.min() != 0 || domain.max() != 1) {
        bez = portion(input, domain) - p;
        partial = true;
    } else {
        bez = input - p;
    }

    // find extrema of the function x(t)^2 + y(t)^2
    // use the fact that (f^2)' = 2 f f'
    // this reduces the order of the distance function by 1
    D2<Bezier> deriv = derivative(bez);
    std::vector<Coord> ts = (multiply(bez[X], deriv[X]) + multiply(bez[Y], deriv[Y])).roots();

    Coord t = -1, mind = infinity();
    for (unsigned i = 0; i < ts.size(); ++i) {
        Coord droot = L2sq(bez.valueAt(ts[i]));
        if (droot < mind) {
            mind = droot;
            t = ts[i];
        }
    }

    // also check endpoints
    Coord dinitial = L2sq(bez.at0());
    Coord dfinal = L2sq(bez.at1());

    if (dinitial < mind) {
        mind = dinitial;
        t = 0;
    }
    if (dfinal < mind) {
        //mind = dfinal;
        t = 1;
    }

    if (partial) {
        t = domain.valueAt(t);
    }
    return t;
}

////////////////////////////////////////////////////////////////////////////////
// D2<SBasis> versions

/*
 * Return the parameter t of the nearest time value on the portion of the curve "c",
 * related to the interval [from, to], to the point "p".
 * The needed curve derivative "dc" is passed as parameter.
 * The function return the first nearest time value to "p" that is found.
 */

double nearest_time(Point const& p,
                    D2<SBasis> const& c,
                    D2<SBasis> const& dc,
                    double from, double to )
{
    if ( from > to ) std::swap(from, to);
    if ( from < 0 || to > 1 )
    {
        THROW_RANGEERROR("[from,to] interval out of bounds");
    }
    if (c.isConstant()) return from;
    SBasis dd = dot(c - p, dc);
    //std::cout << dd << std::endl;
    std::vector<double> zeros = Geom::roots(dd);

    double closest = from;
    double min_dist_sq = L2sq(c(from) - p);
    for ( size_t i = 0; i < zeros.size(); ++i )
    {
        double distsq = L2sq(c(zeros[i]) - p);
        if ( min_dist_sq > L2sq(c(zeros[i]) - p) )
        {
            closest = zeros[i];
            min_dist_sq = distsq;
        }
    }
    if ( min_dist_sq > L2sq( c(to) - p ) )
        closest = to;
    return closest;

}

/*
 * Return the parameters t of all the nearest points on the portion of
 * the curve "c", related to the interval [from, to], to the point "p".
 * The needed curve derivative "dc" is passed as parameter.
 */

std::vector<double>
all_nearest_times(Point const &p,
                  D2<SBasis> const &c,
                  D2<SBasis> const &dc,
                  double from, double to)
{
    if (from > to) {
        std::swap(from, to);
    }
    if (from < 0 || to > 1) {
        THROW_RANGEERROR("[from,to] interval out of bounds");
    }

    std::vector<double> result;
    if (c.isConstant()) {
        result.push_back(from);
        return result;
    }
    SBasis dd = dot(c - p, dc);

    std::vector<double> zeros = Geom::roots(dd);
    std::vector<double> candidates;
    candidates.push_back(from);
    candidates.insert(candidates.end(), zeros.begin(), zeros.end());
    candidates.push_back(to);
    std::vector<double> distsq;
    distsq.reserve(candidates.size());
    for (unsigned i = 0; i < candidates.size(); ++i) {
        distsq.push_back(L2sq(c(candidates[i]) - p));
    }
    unsigned closest = 0;
    double dsq = distsq[0];
    for (unsigned i = 1; i < candidates.size(); ++i) {
        if (dsq > distsq[i]) {
            closest = i;
            dsq = distsq[i];
        }
    }
    for (unsigned i = 0; i < candidates.size(); ++i) {
        if (distsq[closest] == distsq[i]) {
            result.push_back(candidates[i]);
        }
    }
    return result;
}


////////////////////////////////////////////////////////////////////////////////
// Piecewise< D2<SBasis> > versions


double nearest_time(Point const &p,
                    Piecewise< D2<SBasis> > const &c,
                    double from, double to)
{
    if (from > to) std::swap(from, to);
    if (from < c.cuts[0] || to > c.cuts[c.size()]) {
        THROW_RANGEERROR("[from,to] interval out of bounds");
    }

    unsigned si = c.segN(from);
    unsigned ei = c.segN(to);
    if (si == ei) {
        double nearest =
            nearest_time(p, c[si], c.segT(from, si), c.segT(to, si));
        return c.mapToDomain(nearest, si);
    }

    double t;
    double nearest = nearest_time(p, c[si], c.segT(from, si));
    unsigned int ni = si;
    double dsq;
    double mindistsq = distanceSq(p, c[si](nearest));
    Rect bb;
    for (unsigned i = si + 1; i < ei; ++i) {
        bb = *bounds_fast(c[i]);
        dsq = distanceSq(p, bb);
        if ( mindistsq <= dsq ) continue;

        t = nearest_time(p, c[i]);
        dsq = distanceSq(p, c[i](t));
        if (mindistsq > dsq) {
            nearest = t;
            ni = i;
            mindistsq = dsq;
        }
    }
    bb = *bounds_fast(c[ei]);
    dsq = distanceSq(p, bb);
    if (mindistsq > dsq) {
        t = nearest_time(p, c[ei], 0, c.segT(to, ei));
        dsq = distanceSq(p, c[ei](t));
        if (mindistsq > dsq) {
            nearest = t;
            ni = ei;
        }
    }
    return c.mapToDomain(nearest, ni);
}



std::vector<double>
all_nearest_times(Point const &p,
                  Piecewise< D2<SBasis> > const &c,
                  double from, double to)
{
    if (from > to) {
        std::swap(from, to);
    }
    if (from < c.cuts[0] || to > c.cuts[c.size()]) {
        THROW_RANGEERROR("[from,to] interval out of bounds");
    }

    unsigned si = c.segN(from);
    unsigned ei = c.segN(to);
    if ( si == ei )
    {
        std::vector<double>	all_nearest =
            all_nearest_times(p, c[si], c.segT(from, si), c.segT(to, si));
        for ( unsigned int i = 0; i < all_nearest.size(); ++i )
        {
            all_nearest[i] = c.mapToDomain(all_nearest[i], si);
        }
        return all_nearest;
    }
    std::vector<double> all_t;
    std::vector< std::vector<double> > all_np;
    all_np.push_back( all_nearest_times(p, c[si], c.segT(from, si)) );
    std::vector<unsigned> ni;
    ni.push_back(si);
    double dsq;
    double mindistsq = distanceSq( p, c[si](all_np.front().front()) );
    Rect bb;

    for (unsigned i = si + 1; i < ei; ++i) {
        bb = *bounds_fast(c[i]);
        dsq = distanceSq(p, bb);
        if ( mindistsq < dsq ) continue;
        all_t = all_nearest_times(p, c[i]);
        dsq = distanceSq( p, c[i](all_t.front()) );
        if ( mindistsq > dsq )
        {
            all_np.clear();
            all_np.push_back(all_t);
            ni.clear();
            ni.push_back(i);
            mindistsq = dsq;
        }
        else if ( mindistsq == dsq )
        {
            all_np.push_back(all_t);
            ni.push_back(i);
        }
    }
    bb = *bounds_fast(c[ei]);
    dsq = distanceSq(p, bb);
    if (mindistsq >= dsq) {
        all_t = all_nearest_times(p, c[ei], 0, c.segT(to, ei));
        dsq = distanceSq( p, c[ei](all_t.front()) );
        if (mindistsq > dsq) {
            for (unsigned int i = 0; i < all_t.size(); ++i) {
                all_t[i] = c.mapToDomain(all_t[i], ei);
            }
            return all_t;
        } else if (mindistsq == dsq) {
            all_np.push_back(all_t);
            ni.push_back(ei);
        }
    }
    std::vector<double> all_nearest;
    for (unsigned i = 0; i < all_np.size(); ++i) {
        for (unsigned int j = 0; j < all_np[i].size(); ++j) {
            all_nearest.push_back( c.mapToDomain(all_np[i][j], ni[i]) );
        }
    }
    all_nearest.erase(std::unique(all_nearest.begin(), all_nearest.end()),
                      all_nearest.end());
    return all_nearest;
}

namespace {
double clip(double t, double lower = 0, double upper = 1) {
    if (t < lower) {
        return lower;
    }
    if (t > upper) {
        return upper;
    }
    return t;
}

/**
 * @brief nearest_times Calculates the pair of times where two curves are closest to each other
 * @param [inout] result
 * @param stepsize Initial guess for the stepsize, is overwritten if <= 0.
 * @param dist Calculated distance between the curves (optional)
 * @param max_it Maximum number of iterations used.
 */
void nearest_times(std::pair<Coord, Coord>& result,
                   Curve const& a,
                   Curve const& b,
                   double& stepsize,
                   Coord * dist = NULL,
                   size_t const max_it = 20) {
    Coord& t_a = result.first;
    Coord& t_b = result.second;

    Curve const* d_a = a.derivative();
    Curve const* d_b = b.derivative();

    Point residual = a.pointAt(t_a) - b.pointAt(t_b);
    double best_distance = dot(residual, residual);

    if (stepsize <= 0) {
        stepsize = 1.0 / (std::max(
                                     std::abs(dot(residual, d_a->pointAt(t_a))),
                                     std::abs(dot(residual, d_b->pointAt(t_b)))
                                     ));
    }

    for (size_t ii = 0; ii < max_it; ++ii) {
        residual = a.pointAt(t_a) - b.pointAt(t_b);
        double const gradient_a =  dot(residual, d_a->pointAt(t_a));
        double const gradient_b = -dot(residual, d_b->pointAt(t_b));

        double const new_a = clip(t_a - stepsize * gradient_a);
        double const new_b = clip(t_b - stepsize * gradient_b);
        double const update_a = new_a - t_a;
        double const update_b = new_b - t_b;

        double const new_distance = distanceSq(a.pointAt(new_a), b.pointAt(new_b));
        if (new_distance < best_distance) {
            best_distance = new_distance;
            t_a = new_a;
            t_b = new_b;
        }
        else {
            stepsize /= 2;
        }
        if (update_a * update_a + update_b * update_b < 1e-32) {
            break;
        }
    }

    if (dist) {
        *dist = std::sqrt(best_distance);
    }
}

} // anonymous namespace

std::pair<Coord, Coord> nearest_times(
        Curve const& a,
        Curve const& b,
        Coord * dist) {

    std::pair<Coord, Coord> result(0,0);
    double best_distance = distance(a.pointAt(0), b.pointAt(0));
    double best_stepsize = -1;
    std::vector<std::pair<Coord, Coord> > initial_guesses;
    size_t ig_grid_size = 4;
    for (size_t ii = 0; ii <= ig_grid_size; ++ii) {
        double const t_a = (double)ii / ig_grid_size;
        for (size_t jj = 0; jj <= ig_grid_size; ++jj) {
            double const t_b = (double) jj / ig_grid_size;
            initial_guesses.push_back(std::pair<Coord, Coord>(t_a, t_b));
        }
    }

    for (std::pair<Coord, Coord>& ig : initial_guesses) {
        double current_distance = std::numeric_limits<Coord>::max();
        double current_stepsize = -1;
        nearest_times(ig, a, b, current_stepsize, &current_distance);
        if (current_distance < best_distance) {
            best_distance = current_distance;
            best_stepsize = current_stepsize;
            result = ig;
        }
    }

    nearest_times(result, a, b, best_stepsize, &best_distance, 1000);


    if (dist) {
        *dist = best_distance;
    }
    return result;
}

std::pair<PathTime, PathTime> nearest_times(
        Path const& a,
        Path const& b,
        Coord *dist) {
    double min_dist = std::numeric_limits<double>::max();
    std::pair<PathTime, PathTime> result;
    PathTime& time_a = result.first;
    PathTime& time_b = result.second;

    for (size_t ii = 0; ii < a.size(); ++ii) {
        Curve const& curve_a = a[ii];
        Rect const bounds_a = curve_a.boundsFast();
        for (size_t jj = 0; jj < b.size(); ++jj) {
            Curve const& curve_b = b[jj];
            if (distance(bounds_a, curve_b.boundsFast()) > min_dist) {
                continue;
            }
            double current_dist = std::numeric_limits<double>::max();
            std::pair<Coord, Coord> current_result = nearest_times(curve_a, curve_b, &current_dist);
            if (current_dist < min_dist) {
                min_dist = current_dist;
                time_a.curve_index = ii;
                time_a.t = current_result.first;
                time_b.curve_index = jj;
                time_b.t = current_result.second;
            }
        }
    }
    if (dist) {
        *dist = min_dist;
    }
    return result;
}

} // end namespace Geom


