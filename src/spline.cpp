#include "Spline.hpp"

int main(int argc, char *argv[]) {
    QuadraticSplineSegment segment;

    segment.setCtrlPointsAbscissa(0, 1);
    segment.setCtrlPointsOrdinate(0, 0);
    segment.setTangent(0);
    segment.computeParams();
    std::cout << segment.evalFunction(0) << std::endl;
    std::cout << segment.evalDerivative(0) << std::endl;
    std::cout << segment.evalFunction(1) << std::endl;
    std::cout << segment.evalDerivative(1) << std::endl;
    std::cout << segment.arcLength() << std::endl;
    return 0;
}