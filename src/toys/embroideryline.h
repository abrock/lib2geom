#ifndef EMBROIDERYLINE_H
#define EMBROIDERYLINE_H

#include <2geom/d2.h>
#include <2geom/intersection-graph.h>
#include <2geom/path.h>
#include <2geom/sbasis.h>
#include <2geom/svg-path-parser.h>
#include <2geom/transforms.h>
#include <2geom/sbasis-to-bezier.h>
#include <2geom/svg-path-writer.h>
#include <2geom/curve.h>

#include <vector>

using namespace Geom;

class EmbroideryLine;

class OutlineIntersection : public Point {
public:
    /**
     * @brief time The Pathtime of the intersection wrt. the outline path.
     */
    PathTime time;

    /**
     * @brief height The height of the intersection, i.e. the index of the original curve
     */
    size_t height = -1;

    /**
     * @brief outlineStitchHeight The index of the corresponding entry in the outlineStitches vector.
     */
    size_t outline_stitch_index = -1;

    /**
     * @brief index Index of the intersection in the ordered set of intersections with the outline.
     * This allows to check if two intersections are next to each other
     */
    size_t index = 0;

    bool optional = true;

    EmbroideryLine* line;

    OutlineIntersection(
            const Point p = Point()
            , const PathTime _time = PathTime()
            , const size_t _height = -1
            , const size_t _outline_stitch_index = -1
            , const bool _optional = true)
        : Point(p)
        , time(_time)
        , height(_height)
        , outline_stitch_index(_outline_stitch_index)
        , optional(_optional) {}

    void setPoint(Point p) {
        x() = p.x();
        y() = p.y();
    }

    void setPoint(Path path, PathTime time) {
        setPoint(path.pointAt(time));
    }

    friend bool operator< (const OutlineIntersection &c1, const OutlineIntersection &c2) {
        return c1.time < c2.time;
    }

    friend bool operator == (const OutlineIntersection &c1, const OutlineIntersection &c2) {
        return c1.time == c2.time;
    }

    friend bool operator != (const OutlineIntersection &c1, const OutlineIntersection &c2) {
        return !(c1.time == c2.time);
    }

    friend std::ostream& operator << (std::ostream& out, const OutlineIntersection& obj) {
        out << "Time: " << obj.time << ", height: " << obj.height;
        return out;
    }

};

class EmbroideryLine : public std::vector<Point> {
public:
    /**
     * @brief index Number of the curves this EmbroideryLine belongs to.
     */
    size_t level;

    /**
     * @brief up The index of the closest EmbroideryLine in the next curve.
     */
    size_t up;

    /**
     * @brief down All indices of EmbroideryLines on the previous curves.
     */
    std::vector<size_t> down;

    /**
     * @brief startInter PathTime of the point on the outline where this line starts.
     */
    OutlineIntersection startInter;

    /**
     * @brief endInter PathTime of the point on the outline where this line ends.
     */
    OutlineIntersection endInter;

    EmbroideryLine(
            const std::vector<Point>& stitches,
            const size_t _level,
            const OutlineIntersection& _startInter,
            const OutlineIntersection& _endInter) :

        level(_level),
        startInter(_startInter),
        endInter(_endInter) {
        insert(begin(), stitches.begin(), stitches.end());
        //std::cout << "Constructing EmbroideryLine with start: " << startInter << ", end: " << endInter << std::endl;
    }

    friend std::ostream& operator << (std::ostream& out, const EmbroideryLine& obj) {
        out << " Level: " << obj.level
            << " Start: " << obj.startInter
            << " Stop:  " << obj.endInter << std::endl;
        return out;
    }


};

class EmbroideryLineLevel : public std::vector<EmbroideryLine> {
public:
    size_t height;
};


#endif // EMBROIDERYLINE_H
