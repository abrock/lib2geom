/*
 * Circle and Elliptical Arc Fitting Example
 *
 * Authors:
 *      Marco Cecchetti <mrcekets at gmail.com>
 *
 * Copyright 2008  authors
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

#include <memory>
#include <2geom/circle.h>
#include <2geom/elliptical-arc.h>
#include <2geom/numeric/fitting-tool.h>
#include <2geom/numeric/fitting-model.h>

#include <toys/path-cairo.h>
#include <toys/toy-framework-2.h>


using namespace Geom;

class PathPathNearest : public Toy
{
private:
    void draw( cairo_t *cr, std::ostringstream *notify,
               int width, int height, bool save, std::ostringstream *timer_stream)
    {
        if (first_time)
        {
            first_time = false;
        }

        cairo_set_source_rgba(cr, 0.3, 0.3, 0.3, 1.0);
        cairo_set_line_width (cr, 0.3);

        CubicBezier a(psh.pts[0], psh.pts[1], psh.pts[2], psh.pts[3]);
        Path a_path;
        a_path.append(a);

        CubicBezier b(psh.pts[4], psh.pts[5], psh.pts[6], psh.pts[7]);
        Path b_path;
        b_path.append(b);

        cairo_path(cr, a_path);
        cairo_path(cr, b_path);



        Timer tm;
        tm.start();
        double min_dist = 0;
        std::pair<Coord, Coord> times = nearest_times(a_path[0], b_path[0], &min_dist);
        long long time = 0;
        tm.lap(time);
        *notify << time/1000 << "us";
        Point point_a = a_path.pointAt(times.first);
        Point point_b = b_path.pointAt(times.second);

        Path line;
        line.append(LineSegment(point_a, point_b));
        cairo_path(cr, line);

        Rect rect_a(psh.pts[8] - Geom::Point(40,60), psh.pts[8] + Geom::Point(40,60));
        cairo_rectangle(cr, rect_a);

        Rect rect_b(psh.pts[9] - Geom::Point(60,20), psh.pts[9] + Geom::Point(60,20));
        cairo_rectangle(cr, rect_b);

        *notify << std::endl << "Rect dist: " << distance(rect_a, rect_b);
        *notify << std::endl << "Curve dist: " << min_dist;


        cairo_stroke(cr);

        Toy::draw(cr, notify, width, height, save,timer_stream);
    }

public:
    PathPathNearest()
    {
        first_time = true;

        psh.pts.resize(10);
        // First bezier curve
        psh.pts[0] = Point(250, 250);
        psh.pts[1] = Point(300, 200);
        psh.pts[2] = Point(350, 300);
        psh.pts[3] = Point(400, 250);

        // Second bezier curve
        psh.pts[4] = Point(250, 400);
        psh.pts[5] = Point(325, 450);
        psh.pts[6] = Point(325, 350);
        psh.pts[7] = Point(400, 400);

        // middle points of two rectangles for testing distance(Rect, Rect)
        psh.pts[8] = Point(600, 400);
        psh.pts[9] = Point(600, 500);

        for (Point& p : psh.pts) {
            p -= Point(150,150);
            p *= 2;
        }

        handles.push_back(&psh);
    }

private:
    bool first_time;
    PointSetHandle psh;
};



int main(int argc, char **argv)
{
    init( argc, argv, new PathPathNearest(), 600, 600 );
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
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
