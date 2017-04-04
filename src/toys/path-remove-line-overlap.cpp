/*
 * Authors:
 *   Alexander Brock <Brock.Alexander@web.de>
 *
 * Copyright 2017 Authors
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
#include <2geom/d2.h>
#include <2geom/sbasis.h>
#include <2geom/bezier-to-sbasis.h>
#include <2geom/bezier.h>
#include <2geom/sbasis-geometric.h>

#include <toys/path-cairo.h>
#include <toys/toy-framework-2.h>
#include "path-cairo.h"
#include <cairo/cairo.h>
#include <2geom/svg-path-writer.h>
#include <2geom/bezier-utils.h>

#include <vector>
using std::vector;
using namespace Geom;

double asymmetric_furthest_distance(Path const& a, Path const& b) {

    double result = 0;
    for (size_t ii = 0; ii < a.size(); ++ii) {
        Curve const& curve = a[ii];
        for (double t = 0; t <= 1; t += .01) {
            double current_dist = 0;
            b.nearestTime(curve.pointAt(t), &current_dist);
            result = std::max(result, current_dist);
        }
    }

    return result;
}

double symmetric_furthest_distance(Path const& a, Path const& b) {
    return std::max(
                asymmetric_furthest_distance(a,b),
                asymmetric_furthest_distance(b,a)
                );
}

void plot_bezier_with_handles(cairo_t *cr, BezierCurve const& bez, Geom::Point offset = Geom::Point(0,0), double size = 10) {
    Path tmp_path;
    Affine translation;
    translation.setTranslation(offset);
    tmp_path.append(bez);
    tmp_path *= translation;

    cairo_path(cr, tmp_path);
    draw_line_seg(cr, offset + bez[0], offset + bez[1]);
    draw_line_seg(cr, offset + bez[bez.size()-2], offset + bez[bez.size()-1]);
    for (size_t ii = 0; ii < bez.size(); ++ii) {
        Geom::Point const dx(size/2, 0);
        Geom::Point const dy(0, size/2);
        Geom::Point point = bez[ii];
        draw_line_seg(cr, offset + point + dx, offset + point - dx);
        draw_line_seg(cr, offset + point + dy, offset + point - dy);
    }
}

void plot_bezier_with_handles(cairo_t *cr, Path const& path, Geom::Point offset = Geom::Point(0,0), double size = 10) {
    for (size_t ii = 0; ii < path.size(); ++ii) {
        Curve const& curve = path[ii];
        BezierCurve const& bez = dynamic_cast<BezierCurve const&>(curve);
        plot_bezier_with_handles(cr, bez, offset, size);
    }
}

Geom::Path subdivide(CubicBezier const& bez, std::vector<double> const& times_in) {
    Path result;
    // First we need to sort the times ascending.
    std::vector<CubicBezier> curves = bez.subdivide(times_in);
    for (size_t ii = 0; ii < curves.size(); ++ii) {
        result.append(curves[ii]);
    }
    return result;
}

#define SIZE 14

class BezierFitTester: public Toy {
public:
    PointSetHandle b_handle;
    void draw(cairo_t *cr,
              std::ostringstream *notify,
              int width, int height, bool save, std::ostringstream *timer_stream) {

        if (first_time)
        {
            first_time = false;
        }

        Path original_path;
        original_path.append(CubicBezier(b_handle.pts[0], b_handle.pts[1], b_handle.pts[2], b_handle.pts[3]));
        original_path.append(QuadraticBezier(b_handle.pts[3], b_handle.pts[4], b_handle.pts[5]));
        original_path.append(LineSegment(b_handle.pts[5], b_handle.pts[6]));

        Path overlap_path;
        overlap_path.append(CubicBezier(b_handle.pts[7], b_handle.pts[8], b_handle.pts[9], b_handle.pts[10]));
        overlap_path.append(QuadraticBezier(b_handle.pts[10], b_handle.pts[11], b_handle.pts[12]));
        overlap_path.append(LineSegment(b_handle.pts[12], b_handle.pts[13]));
        overlap_path.close();

        cairo_set_line_width (cr, 2.);
        cairo_stroke(cr);

        setRainbowColor(cr, 0);
        plot_bezier_with_handles(cr, original_path);

        setRainbowColor(cr, 1);
        plot_bezier_with_handles(cr, overlap_path);

        {
            PathVector pieces = PathVector(original_path).removeLineOverlap(overlap_path);
            size_t counter = 0;
            for (PathVector::const_iterator it = pieces.begin(); it != pieces.end(); ++it) {
                setRainbowColor(cr, counter++);
                plot_bezier_with_handles(cr, *it, Point(0, 350));
            }
        }


        cairo_stroke(cr);

        Toy::draw(cr, notify, width, height, save,timer_stream);
    }

public:
    BezierFitTester(){
        first_time = true;
        for(int i = 0; i < SIZE; i++) {
            b_handle.push_back(150+uniform()*300,150+uniform()*300);
        }
        size_t index = 0;
        b_handle.pts[index++] = Geom::Point(70,250);
        b_handle.pts[index++] = Geom::Point(200,150);
        b_handle.pts[index++] = Geom::Point(200,350);
        b_handle.pts[index++] = Geom::Point(350,200);
        b_handle.pts[index++] = Geom::Point(450,100);
        b_handle.pts[index++] = Geom::Point(550,200);
        b_handle.pts[index++] = Geom::Point(650,200);

        b_handle.pts[index++] = Geom::Point(70,350);
        b_handle.pts[index++] = Geom::Point(210,160);
        b_handle.pts[index++] = Geom::Point(230,160);
        b_handle.pts[index++] = Geom::Point(350,300);
        b_handle.pts[index++] = Geom::Point(450,200);
        b_handle.pts[index++] = Geom::Point(550,300);
        b_handle.pts[index++] = Geom::Point(650,300);
        handles.push_back(&b_handle);
    }
private:
    std::vector<Slider> sliders;
    bool first_time;

    void setRainbowColor(cairo_t* cr, size_t index) {
        switch(index) {
        case 0:
            cairo_set_source_rgba (cr, 100./255, 39./255, 134./255, 1);
            break;
        case 1:
            cairo_set_source_rgba (cr, 20./255, 151./255, 214./255, 1);
            break;
        case 2:
            cairo_set_source_rgba (cr, 92./255, 187./255, 70./255, 1);
            break;
        case 3:
            cairo_set_source_rgba (cr, 244./255, 108./255, 47./255, 1);
            break;
        case 4:
            cairo_set_source_rgba (cr, 239./255, 60./255, 63./255, 1);
            break;
        default:
            cairo_set_source_rgba (cr, 105./255, 171./255, 222./255, 1);
        }
    }
};

int main(int argc, char **argv) {
    init(argc, argv, new BezierFitTester);
    return 0;
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
//vim:filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99:
