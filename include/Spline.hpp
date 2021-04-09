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
        std::cout << a2 << " " << a1 << " " << a0 << std::endl;
    }

    void getParams(double &a_2, double &a_1, double &a_0) const {
        a_2 = a2;
        a_1 = a1;
        a_0 = a0;
    }

    double getTangent() const {
        return t1;
    }

    double arcLength(double gamma) {
        computeParams();
        double yi = y1, yi1, xi = x1, xi1, len, ti;
        while (norm(xi - x2, yi - y2) > gamma * 3) {
            ti = evalDerivative(xi);
            yi = evalFunction(xi);
            xi1 = xi + gamma / norm(1, ti);
            yi1 = evalFunction(xi1);
            len += norm(xi1 - xi, yi1 - yi);
            std::cout << ((norm(xi - x2, yi - y2) - gamma) / gamma - 1.0) * 100 << " " << xi << std::endl;
            xi = xi1;
        }
        return len;
    }

    double evalFunction(double x) const {
        return a2 * pow(x, 2) + a1 * x + a0;
    }

    double evalDerivative(double x) const {
        return a2 * x + a1;
    }

protected:

    static inline double norm(double x, double y) {
        return sqrt(pow(x, 2) + pow(y, 2));
    }

    double t1, x1, x2, y1, y2;
    double a2, a1, a0;
};

class QuadraticSplineCurve {
public:

    explicit QuadraticSplineCurve(std::vector<std::vector<double>> coordinates, double tangent) : coord(
            std::move(coordinates)), t(tangent) {}

    void compute() {
        QuadraticSplineSegment segment{};
        double x1, y1, x2, y2;

        for (int i = 0; i < coord.size() - 1; ++i) {
            auto pts1 = coord[i];
            auto pts2 = coord[i + 1];
            x1 = pts1[0];
            y1 = pts1[1];
            x2 = pts2[0];
            y2 = pts2[1];

            segment.setCtrlPointsAbscissa(x1, x2);
            segment.setCtrlPointsOrdinate(y1, y2);

            if (i == 0) {
                segment.setTangent(t);
            } else {
                segment.setTangent(curve[curve.size() - 1].evalDerivative(x2));
            }

            segment.computeParams();

            curve.push_back(segment);
        }
    }

    double evalFunction(double x) {
        double x1, x2;

        for (int i = 0; i < coord.size() - 1; ++i) {

            auto pts1 = coord[i];
            auto pts2 = coord[i + 1];

            x1 = pts1[0];
            x2 = pts2[0];

            if (x1 < x2) {
                if ((x >= x1) and (x <= x2)) {
                    return curve[i].evalFunction(x);
                }
            } else {
                if ((x <= x1) and (x >= x2)) {
                    return curve[i].evalFunction(x);
                }
            }
        }
    }

    double length(double gamma) {
        cLength = 0;
        for (auto segment : curve) {
            cLength += segment.arcLength(gamma/1000.0);
        }
    }

protected:
    std::vector<std::vector<double>> coord;
    std::vector<QuadraticSplineSegment> curve;
    double t, cLength{};
};

#endif //SPLINE_SPLINE_HPP
