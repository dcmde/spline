#ifndef SPLINE_SPLINE_HPP
#define SPLINE_SPLINE_HPP

#include <Eigen/Dense>
#include <utility>
#include <vector>
#include <iostream>

class QuadraticSplineSegment {
public:

    void setCtrlPointsAbscissa(double x_1, double x_2) {
        x1 = x_1;
        x2 = x_2;
    }

    void setCtrlPointsOrdinate(double y_1, double y_2) {
        y1 = y_1;
        y2 = y_2;
    }

    void setTangent(double tangent) {
        t1 = tangent;
    }

    bool computeParams() {
        Eigen::Matrix3d A;
        Eigen::Vector3d Y, X;

        A << pow(x1, 2), x1, 1,
                pow(x2, 2), x2, 1,
                2 * x1, 1, 0;

        Y << y1, y2, t1;

        X = A.colPivHouseholderQr().solve(Y);

        a2 = X[0];
        a1 = X[1];
        a0 = X[2];
    }

    void getParams(double &a_2, double &a_1, double &a_0) const {
        a_2 = a2;
        a_1 = a1;
        a_0 = a0;
    }

    double arcLength() const {
        double temp1, temp2, temp3;
        if (a2 != 0) {
            temp2 = sqrt(pow((2 * a2 * x2 + a1 + 1), 3)) / (3 * a2);
            temp1 = sqrt(pow((2 * a2 * x1 + a1 + 1), 3)) / (3 * a2);
            temp3 = temp2 - temp1;
        } else {
            temp3 = sqrt(1 + pow(a1, 2)) * (x2 - x1);
        }
        return temp3;
    }

    double evalFunction(double x) const {
        return a2 * pow(x, 2) + a1 * x + a0;
    }

    double evalDerivative(double x) const {
        return a2 * x + a1;
    }

protected:
    double t1, x1, x2, y1, y2;
    double a2, a1, a0;
};

class CubicSplineCurve {
public:
    explicit CubicSplineCurve(std::vector<std::vector<double>> coordinates) : coord(std::move(coordinates)) {}

protected:
    std::vector<std::vector<double>> coord;
};

#endif //SPLINE_SPLINE_HPP
