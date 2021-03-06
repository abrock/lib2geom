#pragma once
#ifndef HAIR_H
#define HAIR_H

#include <2geom/d2.h>
#include <2geom/intersection-graph.h>
#include <2geom/path.h>
#include <2geom/sbasis.h>

#include <2geom/svg-path-parser.h>
#include <2geom/transforms.h>
#include <2geom/sbasis-to-bezier.h>
#include <2geom/svg-path-writer.h>
#include <2geom/curve.h>

//#include <toys/path-cairo.h>
//#include <toys/toy-framework-2.h>

//#include <cairo/cairo.h>
//#include <gtk-2.0/gdk/gdkevents.h>

#include <cstdlib>
#include <boost/algorithm/minmax_element.hpp>

#include <iostream>
#include <fstream>
#include <random>

#include "helper/geom-pathstroke.h"

#include "RunningStats.h"

#include "memoryConsumption.h"

#include "embroideryline.h"

#include "embroideryarea.h"


using namespace Geom;


/*
class OutlineStitch : public Geom::Point {
public:
    PathTime time;
    bool optional;
    OutlineStitch(Geom::Point p, PathTime _time, bool _optional)
        : Geom::Point(p)
        , time(_time)
        , optional(_optional)
    {

    }
    friend bool operator< (const OutlineStitch &c1, const OutlineStitch &c2) {
        return c1.time < c2.time;
    }

    friend bool operator == (const OutlineStitch &c1, const OutlineStitch &c2) {
        return c1.time == c2.time;
    }

    friend bool operator != (const OutlineStitch &c1, const OutlineStitch &c2) {
        return !(c1.time == c2.time);
    }
};
*/

class Hair {
    typedef char BOOL;
private:
    /**
     * @brief outline is the (closed) path which later contains all the hair curves
     */
    Path _outline;
    /**
     * @brief curve The inital curve from which all hair curves are derivated.
     */
    Path _curve;

    Path _feather_template;

    Path _feather_curve;

    /**
     * @brief curves A vector of curves generated by the method. On this curve lies
     */
    PathVector _curves;

    /**
     * @brief offset Relative offset between the stitches of one line and the next. Relative to current stitch length.
     */
    Coord _offset = .311;

    double _line_spacing = .35;

    /**
     * @brief areaGrow Make the area larger.
     */
    double _area_grow = 0.0;

    /**
     * @brief miter Miter limit for the Inkscape::half_outline method. Paths are expected to be smooth so this limit shouldn't do anything at all.
     */
    Coord _miter = 0;

    //Inkscape::LineJoinType _join_type = Inkscape::JOIN_MITER; // 4625
    Inkscape::LineJoinType _join_type = Inkscape::JOIN_BEVEL; // 4625
    // Inkscape::LineJoinType _join_type = Inkscape::JOIN_MITER_CLIP; // 1181
    // Inkscape::LineJoinType _join_type = Inkscape::JOIN_EXTRAPOLATE; // too long

    /**
     * @brief stitchLength Initial stitch length used for the line in the middle.
     * Stitch lenghtes for the other lines are automatically adjusted by the projection
     * of the stitches from one curve to the next.
     */
    Coord _stitch_length = 3.6;

    size_t _max_iter = 600;

    /**
     * @brief stitches Stitches generated by the method.
     */
    std::vector<std::vector<Point> > stitches;

    double _stitch_point_radius = _line_spacing *.45;

    double _line_width = _line_spacing * .5;

    double _max_stitch_length = 10;

    double _min_stitch_length = .1;

    std::vector<double> _initial_stitch_lengths;
    double _initial_stitch_lengths_sum;

    std::mt19937_64 generator;

    double _discrete_outline_param = 2;
    std::vector<Point> _discrete_outline;

    std::vector<std::vector<Point> > _bridges;

    std::vector<EmbroideryLineLevel> _levels;

    std::vector<OutlineIntersection> _intersections;

    std::vector<OutlineIntersection> _outline_stitches;

    std::vector<EmbroideryArea> _areas;

    RunningStats _stitch_length_stats;

    size_t _feather_corner;
    Point _feather_corner_point;

public:

    void setLineSpacing(double const new_spacing);

    void setDensity(double const new_density);


    double meanDist(Path const& template_path, Path const& offsetted, Path const& area_of_interest = Path());

    bool offsetTemplateOnCurve(
            Path const& template_path,
            Path const& last_path,
            Path const& curve,
            PathTime const last_time,
            PathTime & new_time,
            double const target_distance,
            double &offset,
            Path const& area_of_interest = Path());

    void setOutline(Path const& p);

    void setCurve(Path const& p);

    void getFeatherCurves();

    void getFeatherStitches();

    std::vector<Geom::Point> getSingleFeatherStitches(
            Path const& path,
            double const length,
            double offset = 0);

    void setFeatherCurve(Path const& p);

    void setFeatherTemplate(Path const& p);

    void printStats();

    std::string getStats();

    std::vector<Point> greedySolution;

    static Point getCenter(const std::vector<Point>& stitches);

    static Point getCenter(const std::vector<std::vector<Point> >& stitches);

    void assembleGreedySolution();

    void getOutlineStitches();

    void writeOutlineIntersections(std::ostream& out, std::vector<OutlineIntersection> &_intersections);

    static double curveLength(const Geom::Path& _curve,
                              const Geom::PathTime start,
                              const Geom::PathTime stop,
                              const double dt = 1.0/20000);

    /**
     * @brief writeStitches Write the stitches into a text file suitable for libembroidery-convert
     * @param filename
     */
    void writeStitches(const char* filename);

    void writeStitches(const std::string filename);
    /**
     * @brief writeStitches Write the stitches into a text file suitable for libembroidery-convert
     * @param filename
     */
    static void writeStitches(const std::vector<Point>& stitches, const char* filename);

    /**
     * @brief writeStitches Write the stitches into a text file suitable for libembroidery-convert
     * @param filename
     */
    static void writeStitches(const std::vector<std::vector<Point> >& stitches, const char* filename);

    /**
     * @brief writeStitches Write the stitches into a text file suitable for libembroidery-convert
     * @param filename
     */
    static void writeStitches(const std::vector<std::vector<Point> >& stitches, std::string filename);

    /**
     * @brief writeStitches Write the stitches into a text file suitable for libembroidery-convert
     * @param out
     */
    static void writeStitches(const std::vector<std::vector<Point> >& stitches, std::ostream& out);

    /**
     * @brief writeStitches Write the stitches into a text file suitable for libembroidery-convert
     * @param filename
     */
    static void writeStitches(const std::vector<Point>& stitches, const std::string filename);

    static std::string outputStitch(Point p);

    void writeStitches(std::ostream& out);

    static void writeStitches(const std::vector<Point>& stitches, std::ostream& out);

    double getStitchLengthes(const std::vector<Point>& stitches, std::vector<double>& lengths);

    Hair(Path _outline, Path _curve);

    Hair() {}

    void run();
    void runFeather();

    size_t countStitches();

    template<class Inner>
    size_t countStitches(const std::vector<Inner>& elements);

    std::vector<Point> purgeSmallStitches(
            const std::vector<Point>& orig,
            const double minLength);

    std::vector<Point> purgeSmallStitchesSub(
            const std::vector<Point>& orig,
            const double minLength);

    std::vector<Point> subDivideLargeStitches(
            const std::vector<Point>& orig,
            const double maxLength);

    void getBoundaryDiscretization();

    void lengthStats(const std::vector<Point>& points);

    std::vector<std::pair<size_t, size_t> > endpointOutlineMatches(
            const std::vector<Point>& _outline,
            const std::vector<std::vector<Point> >& lines);

    size_t bestMatch(const std::vector<Point>& points, const Point& target);

    size_t moduloDist(const size_t _a, const size_t _b, const size_t _total);

    int moduloDistDirection(const size_t _a, const size_t _b, const size_t _total);

    void assembleStitches();

    void assembleAreas();

    void getOutlineIntermediateStitches();

    void assignOutlineIntersections();

    /**
     * @brief isAdjacentLine checks if two lines should belong to the same EmbroideryArea
     * @param a
     * @param b
     * @return true if a and b should belong to the same EmbroideryArea
     */
    bool isAdjacentLine(const EmbroideryLine& a, const EmbroideryLine& b);

    std::vector<Point> assembleStitches(
            const std::vector<std::vector<Point> >& stitches,
            std::vector<std::vector<Point> >& _bridges,
            const size_t initialLine,
            const bool reverseInitial,
            std::vector<BOOL>& shouldCheck);

    int getNearestLine(
            const Point& p,
            const std::vector<BOOL>& visited,
            const std::vector<std::vector<Point> >& lines);

    int getNextLine(
            const int lastLine,
            const Point& p,
            const std::vector<BOOL>& visited,
            const std::vector<std::vector<Point> >& lines);

    /**
     * @brief purgeOutside Removes all stitches outside the given boundary
     * and replaces stitches crossing the boundary by stitches ending / starting at the boundary.
     */
    void purgeOutside();

    double lineLength(const std::vector<Point>& points);


    void getStitches();

    void enlargeCurves();

    void enlargeCurve(Path& _curve);

    void enlargeCurveEnd(Path& _curve);

    /**
     * @brief atLeastOneInside
     * @param points
     * @param boundary
     * @return
     */
    bool atLeastOneInside(const std::vector<Point>& points, const Path& boundary);

    /**
     * @brief vectorOffset moves every point in a vector in the direction given by the difference between itself and the previous point.
     * @param points The original points
     * @param offset The relative offset, usually between -1 and 1
     * @return The moved points.
     */
    std::vector<Point> vectorOffsetOnCurve(
            const std::vector<Point>& points,
            const Path& _curve,
            const double _offset = 0.0);

    bool isInsidePoint(const Point p);

    std::vector<Point> projection(const std::vector<Point> & src, Path& _curve);

    std::vector<Geom::Point> getStitchesOnCurve(
            Path const& curve,
            double const stitchLength,
            PathTime const& start = PathTime());

    Point getOffsettedPointOnCurve(const Path& _curve, PathTime& t, const double targetOffset);

    void makeAreaLarger(Path& _curve, const double _offset);

    std::vector<Point> getShortestConnection(OutlineIntersection a, OutlineIntersection b) const;

    void getCurves();

    void write(const char* filename);

    void write(const std::string filename);

    void write(std::ostream& out);

    void writeCurves(const char* filename);
    void writeCurves(const std::string filename);
    void writeCurves(std::ostream& out);


    void write(std::ostream& out, Path path, std::string color);

    void writeStitchPoints(std::ostream& out, std::string color);

    void writeAreas(const char* filename);

    void writeForwardAreas(const char* filename);
    void writeReverseAreas(const char* filename);

    void writeArea(std::ofstream& out, const EmbroideryArea& area);

    std::string getColor(size_t pos, size_t numColors);

    void addStartStop();

    static void addStartStop(std::vector<Point>& patch);

    template<class Line>
    void writePatches(std::ostream& out, std::vector<Line>& patches);

    template<class C>
    Path getPath(const std::vector<C>& points);

    void writeCircles(
            std::ostream& out,
            const std::vector<Point> & centers,
            double radius,
            std::string color);



    void writeCircle(std::ostream& out, Point center, double radius, std::string color);

};

#endif // HAIR_H
