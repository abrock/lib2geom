#ifndef EMBROIDERYPATCH_H
#define EMBROIDERYPATCH_H
#include <vector>
#include "embroideryline.h"

/**
 * @brief The EmbroideryArea class represents a sub-area of a larger filled embroidery area
 * which can be filled without using connecting sub-stitches.
 */
class EmbroideryArea : public std::vector<EmbroideryLine> {


public:
    const size_t min_level;
    std::vector<Point> forward_stitches;
    std::vector<Point> reverse_stitches;
    OutlineIntersection start;
    OutlineIntersection stop;
    OutlineIntersection r_start;
    OutlineIntersection r_stop;

    EmbroideryArea(EmbroideryLine first_line, const size_t _min_level)
        : min_level(_min_level)
        , start(first_line.startInter)
        , stop(first_line.startInter)
        , r_start(first_line.endInter)
        , r_stop (first_line.endInter)
    {
        push_back(first_line);
    }

    /**
     * @brief finish Create the stitch lists for forward and reverse
     */
    void finish() {
        for (size_t ii = 1; ii < size(); ++ii) {

        }
    }
};

#endif // EMBROIDERYPATCH_H
