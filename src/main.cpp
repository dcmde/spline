#include <iostream>
#include <fstream>
#include <sstream>
#include "BSplines.hpp"

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

    int cpt = 0, numPts;
    double x, y, dx, dy;
    Vec3 temp2;
    std::string line;
    std::vector<double> temp;
    std::vector<Vec3> ctrlPts;
    std::vector<std::vector<double>> ctrlPoints;

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

    std::ofstream ofstream("output.txt");

    if (!ofstream.is_open()) {
        std::cout << "Cannot open the output file." << std::endl;
        return 0;
    }

    std::cout << "Input control point:" << std::endl;
    std::cout << "x y" << std::endl;

    while (std::getline(ifstream, line)) {
        std::cout << line << std::endl;
        split(line, temp);
        ctrlPts.emplace_back(temp[0], temp[1], 0);
        ++cpt;
        temp.clear();
    }

//    std::cout << "Output" << std::endl;
    ofstream << "a x y dx dy" << std::endl;

    BSpline<Vec3, float> spline(cpt, eOPEN_UNIFORM);

    spline.set_ctrl_points(ctrlPts);

    for (int i = 0; i < numPts; ++i) {
        float alpha = i / (float) numPts;
        temp2 = spline.eval_f(alpha);
        x = temp2.x;
        y = temp2.y;


        temp2 = spline.eval_df(alpha);
        dx = temp2.x;
        dy = temp2.y;
//        std::cout << alpha << " " << x << " " << y << " " << dx << " " << dy << std::endl;
        ofstream << alpha << " " << x << " " << y << " " << dx << " " << dy << std::endl;
    }
}