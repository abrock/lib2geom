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

class OutlineIntersection {
public:
    /**
     * @brief outline Reference to the outline Path
     */
    Path& outline;

    /**
     * @brief time The Pathtime of the intersection wrt. the outline path.
     */
    PathTime time;

    /**
     * @brief height The height of the intersection, i.e. the index of the original curve
     */
    size_t height;

    /**
     * @brief outlineStitchHeight The index of the corresponding entry in the outlineStitches vector.
     */
    size_t outline_stitch_index;

    EmbroideryLine* line;

    Point point() {
        return outline.pointAt(time);
    }

    OutlineIntersection(Path& _outline, const PathTime _time, const size_t _height):
        outline(_outline),
        time(_time),
        height(_height) {}

    OutlineIntersection(Path& _outline) : outline(_outline), time(PathTime()), height(0) {
        time.curve_index = 0;
        time.t = 0.0;
    }

    /**
     * @brief index Index of the intersection in the ordered set of intersections with the outline.
     * This allows to check if two intersections are next to each other
     */
    mutable size_t index = 0;

    OutlineIntersection& operator=(OutlineIntersection b) {
        index = b.index;
        height = b.height;
        time = b.time;
        if (outline != b.outline) {
            //throw std::logic_error(std::string("Invalid attempt to assign OutlineIntersection with different outline member to self\nfile: ") + __FILE__ + "\nline: " + std::to_string(__LINE__) + "\n");
        }
        return *this;
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
