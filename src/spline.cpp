#include <fstream>
#include "Spline.hpp"

void split(const std::string &str, std::vector<double> &vec, char delim = ' ') {
    double val;
    std::stringstream ss(str);
    std::string token;

    while (std::getline(ss, token, delim)) {
        val = std::stod(token);
        vec.push_back(val);
    }
}

int main(int argc, char *argv[]) {
    std::vector<std::vector<double>> coordinates;
    std::vector<double> pts;
    std::ofstream ofstream("spline_output.txt");
    std::string line;
    int numPts;
    double x, dx, y;

    if (argc != 3) {
        std::cout << "Enter : file_name number_of_points" << std::endl;
        return 0;
    }

    std::ifstream ifstream(argv[1]);

    if (!ifstream.is_open()) {
        std::cout << "Cannot open the input file." << std::endl;
        return 0;
    }

    numPts = std::stoi(argv[2]);

    std::cout << "Input control point:" << std::endl;
    std::cout << "x y" << std::endl;

    while (std::getline(ifstream, line)) {
        std::cout << line << std::endl;
        split(line, pts);
        coordinates.push_back(pts);
        pts.clear();
    }

    ofstream << "x y" << std::endl;

    QuadraticSplineCurve curve(coordinates, 0);

    curve.compute();
    dx = abs(coordinates[0][0] - coordinates[coordinates.size() - 1][0]) / numPts;
    x = coordinates[0][0];

    for (int i = 0; i < numPts; ++i) {
        ofstream << x << " " << curve.evalFunction(x) << std::endl;
        x += dx;
    }

}