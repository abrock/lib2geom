#ifndef EMBROIDERYOPTIMIZATION_H
#define EMBROIDERYOPTIMIZATION_H

#include <vector>

class EmbroideryAreaChoice {
public:
    size_t areaIndex;
    bool reverse;

    EmbroideryAreaChoice (size_t _areaIndex, bool _reverse) : areaIndex(_areaIndex), reverse(_reverse) {

    }
};

class EmbroideryAreaSequence : public std::vector<EmbroideryAreaChoice> {

};

#endif // EMBROIDERYOPTIMIZATION_H
