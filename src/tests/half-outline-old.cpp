#include <2geom/affine.h>
#include <2geom/transforms.h>
#include <2geom/svg-path-parser.h>
#include <2geom/sbasis-to-bezier.h>

void tangents_old(Geom::Point tang[2], Geom::Curve const& incoming, Geom::Curve const& outgoing)
{
    Geom::Point tang1 = Geom::unitTangentAt(reverse(incoming.toSBasis()), 0.);
    Geom::Point tang2 = outgoing.unitTangentAt(0.);
    tang[0] = tang1, tang[1] = tang2;
}

// Offsetting a line segment is mathematically stable and quick to do
Geom::LineSegment offset_line_old(Geom::LineSegment const& l, double width)
{
    Geom::Point tang1 = Geom::rot90(l.unitTangentAt(0));
    Geom::Point tang2 = Geom::rot90(unitTangentAt(reverse(l.toSBasis()), 0.));

    Geom::Point start = l.initialPoint() + tang1 * width;
    Geom::Point end = l.finalPoint() - tang2 * width;

    return Geom::LineSegment(start, end);
}

void get_cubic_data_old(Geom::CubicBezier const& bez, double time, double& len, double& rad)
{
    // get derivatives
    std::vector<Geom::Point> derivs = bez.pointAndDerivatives(time, 3);

    Geom::Point der1 = derivs[1]; // first deriv (tangent vector)
    Geom::Point der2 = derivs[2]; // second deriv (tangent's tangent)
    double l = Geom::L2(der1); // length

    len = rad = 0;

    // TODO: we might want to consider using Geom::touching_circle to determine the
    // curvature radius here. Less code duplication, but slower

    if (Geom::are_near(l, 0, 1e-4)) {
        l = Geom::L2(der2);
        Geom::Point der3 = derivs.at(3); // try second time
        if (Geom::are_near(l, 0, 1e-4)) {
            l = Geom::L2(der3);
            if (Geom::are_near(l, 0)) {
                return; // this isn't a segment...
            }
            rad = 1e8;
        } else {
            rad = -l * (Geom::dot(der2, der2) / Geom::cross(der2, der3));
        }
    } else {
        rad = -l * (Geom::dot(der1, der1) / Geom::cross(der1, der2));
    }
    len = l;
}

void offset_cubic_old(Geom::Path& p, Geom::CubicBezier const& bez, double width, double tol, size_t levels)
{
    using Geom::X;
    using Geom::Y;

    Geom::Point start_pos = bez.initialPoint();
    Geom::Point end_pos = bez.finalPoint();

    Geom::Point start_normal = Geom::rot90(bez.unitTangentAt(0));
    Geom::Point end_normal = -Geom::rot90(Geom::unitTangentAt(Geom::reverse(bez.toSBasis()), 0.));

    // offset the start and end control points out by the width
    Geom::Point start_new = start_pos + start_normal*width;
    Geom::Point end_new = end_pos + end_normal*width;

    // --------
    double start_rad, end_rad;
    double start_len, end_len; // tangent lengths
    get_cubic_data_old(bez, 0, start_len, start_rad);
    get_cubic_data_old(bez, 1, end_len, end_rad);

    double start_off = 1, end_off = 1;
    // correction of the lengths of the tangent to the offset
    if (!Geom::are_near(start_rad, 0))
        start_off += width / start_rad;
    if (!Geom::are_near(end_rad, 0))
        end_off += width / end_rad;
    start_off *= start_len;
    end_off *= end_len;
    // --------

    Geom::Point mid1_new = start_normal.ccw()*start_off;
    mid1_new = Geom::Point(start_new[X] + mid1_new[X]/3., start_new[Y] + mid1_new[Y]/3.);
    Geom::Point mid2_new = end_normal.ccw()*end_off;
    mid2_new = Geom::Point(end_new[X] - mid2_new[X]/3., end_new[Y] - mid2_new[Y]/3.);

    // create the estimate curve
    Geom::CubicBezier c = Geom::CubicBezier(start_new, mid1_new, mid2_new, end_new);

    // reached maximum recursive depth
    // don't bother with any more correction
    if (levels == 0) {
        p.append(c);
        return;
    }

    // check the tolerance for our estimate to be a parallel curve
    Geom::Point chk = c.pointAt(.5);
    Geom::Point req = bez.pointAt(.5) + Geom::rot90(bez.unitTangentAt(.5))*width; // required accuracy

    Geom::Point const diff = req - chk;
    double const err = Geom::dot(diff, diff);

    if (err < tol) {
        if (Geom::are_near(start_new, p.finalPoint())) {
            p.setFinal(start_new); // if it isn't near, we throw
        }

        // we're good, curve is accurate enough
        try {
          p.append(c);}
        catch (...) {

        }

        return;
    } else {
        // split the curve in two
        std::pair<Geom::CubicBezier, Geom::CubicBezier> s = bez.subdivide(.5);
        offset_cubic_old(p, s.first, width, tol, levels - 1);
        offset_cubic_old(p, s.second, width, tol, levels - 1);
    }
}

void offset_quadratic_old(Geom::Path& p, Geom::QuadraticBezier const& bez, double width, double tol, size_t levels)
{
    // cheat
    // it's faster
    // seriously
    std::vector<Geom::Point> points = bez.controlPoints();
    Geom::Point b1 = points[0] + (2./3) * (points[1] - points[0]);
    Geom::Point b2 = b1 + (1./3) * (points[2] - points[0]);
    Geom::CubicBezier cub = Geom::CubicBezier(points[0], b1, b2, points[2]);
    offset_cubic_old(p, cub, width, tol, levels);
}

void offset_curve_old(Geom::Path& res, Geom::Curve const* current, double width)
{
    double const tolerance = 0.0025;
    size_t levels = 8;

    if (current->isDegenerate()) return; // don't do anything

    // TODO: we can handle SVGEllipticalArc here as well, do that!

    if (Geom::BezierCurve const *b = dynamic_cast<Geom::BezierCurve const*>(current)) {
        size_t order = b->order();
        switch (order) {
        case 1:
            res.append(offset_line_old(static_cast<Geom::LineSegment const&>(*current), width));
            break;
        case 2: {
            Geom::QuadraticBezier const& q = static_cast<Geom::QuadraticBezier const&>(*current);
            offset_quadratic_old(res, q, width, tolerance, levels);
            break;
        }
        case 3: {
            Geom::CubicBezier const& cb = static_cast<Geom::CubicBezier const&>(*current);
            offset_cubic_old(res, cb, width, tolerance, levels);
            break;
        }
        default: {
            Geom::Path sbasis_path = Geom::cubicbezierpath_from_sbasis(current->toSBasis(), tolerance);
            for (size_t i = 0; i < sbasis_path.size(); ++i)
                offset_curve_old(res, &sbasis_path[i], width);
            break;
        }
        }
    } else {
        Geom::Path sbasis_path = Geom::cubicbezierpath_from_sbasis(current->toSBasis(), 0.1);
        for (size_t i = 0; i < sbasis_path.size(); ++i)
            offset_curve_old(res, &sbasis_path[i], width);
    }
}

Geom::Path half_outline_old(Geom::Path const& input, double width, double miter, Inkscape::LineJoinType join = Inkscape::JOIN_BEVEL)
{
    Geom::Path res;
    if (input.size() == 0) return res;

    Geom::Point tang1 = input[0].unitTangentAt(0);
    Geom::Point start = input.initialPoint() + tang1 * width;
    Geom::Path temp;
    Geom::Point tang[2];

    res.setStitching(true);
    temp.setStitching(true);

    res.start(start);

    // Do two curves at a time for efficiency, since the join function needs to know the outgoing curve as well
    const size_t k = (input.back_closed().isDegenerate() && input.closed())
            ?input.size_default()-1:input.size_default();
    for (size_t u = 0; u < k; u += 2) {
        temp.clear();

        offset_curve_old(temp, &input[u], width);

        // on the first run through, there isn't a join
        if (u == 0) {
            res.append(temp);
        } else {
            tangents_old(tang, input[u-1], input[u]);
            outline_join(res, temp, tang[0], tang[1], width, miter, join);
        }

        // odd number of paths
        if (u < k - 1) {
            temp.clear();
            offset_curve_old(temp, &input[u+1], width);
            tangents_old(tang, input[u], input[u+1]);
            outline_join(res, temp, tang[0], tang[1], width, miter, join);
        }
    }

    if (input.closed()) {
        Geom::Curve const &c1 = res.back();
        Geom::Curve const &c2 = res.front();
        temp.clear();
        temp.append(c1);
        Geom::Path temp2;
        temp2.append(c2);
        tangents_old(tang, input.back(), input.front());
        outline_join(temp, temp2, tang[0], tang[1], width, miter, join);
        res.erase(res.begin());
        res.erase_last();
        //
        res.append(temp);
        res.close();
    }

    return res;
}
