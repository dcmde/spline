# BSpline and QuadraticSpline

The BSpline is based on the code of Arkan (https://github.com/brainexcerpts/BSpline). The `bsplines.hpp` allows to compute
bspline functions and the first derivative. The `spline.hpp` allows to compute a quadratic spline function.

## How To

If you want to use the CLI and the jupyter notebook do as follow,

```
mkdir build
cd build
cmake ..
make
```

Then you need to create a file with control points the first column is the `x` coordinates, the second the `y`
coordinates. The separation parameter is a single space. Then you can use the CLI by specifying the file name and the
number of points like this,

```
./spline file.txt 100
./bspline file.txt 100
```

The script will give (`a`, `x`, `y`, `dx`, `dy`) and store the data in the `x_output.txt` file.