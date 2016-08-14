#ifndef EMBROIDERYPATCH_H
#define EMBROIDERYPATCH_H
#include <vector>
#include "embroideryline.h"

/**
 * @brief The EmbroideryArea class represents a sub-area of a larger filled embroidery area
 * which can be filled without using connecting sub-stitches.
 */
class EmbroideryArea : public std::vector<EmbroideryLine> {
    const size_t min_level;
    EmbroideryArea(const size_t _min_level) : min_level(_min_level) {}
};

#endif // EMBROIDERYPATCH_H
