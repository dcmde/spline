# Spline

This repo is based on the code of Arkan. The `splines.hpp` allows to compute spline function.

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
```

Data are stored in the `output.txt` file.