#ifndef SPLINES_HPP
#define SPLINES_HPP

#include <vector>
#include <cmath>


struct Vec3 {

    float x, y, z;

    // -------------------------------------------------------------------------
    /// @name Constructors
    // -------------------------------------------------------------------------

    Vec3() {
        x = float(0);
        y = float(0);
        z = float(0);
    }

    Vec3(float x_, float y_, float z_) {
        x = x_;
        y = y_;
        z = z_;
    }

    explicit Vec3(float v) {
        x = v;
        y = v;
        z = v;
    }

    //Vec3(const Vec3& v) { x = float(v.x); y = float(v.y); z = float(v.z); }

    static Vec3 unit_x() { return Vec3(float(1), float(0), float(0)); }

    static Vec3 unit_y() { return Vec3(float(0), float(1), float(0)); }

    static Vec3 unit_z() { return Vec3(float(0), float(0), float(1)); }

    static Vec3 zero() { return Vec3(float(0), float(0), float(0)); }

    static Vec3 unit_scale() { return Vec3(float(1), float(1), float(1)); }


    void set(float x_, float y_, float z_) {
        x = x_;
        y = y_;
        z = z_;
    }

    // -------------------------------------------------------------------------
    /// @name Overload operators
    // -------------------------------------------------------------------------

    // ----------
    // Additions
    // ----------


    Vec3 operator+(const Vec3 &v_) const {
        return Vec3(x + v_.x, y + v_.y, z + v_.z);
    }

    Vec3 &operator+=(const Vec3 &v_) {
        x += v_.x;
        y += v_.y;
        z += v_.z;
        return *this;
    }

    Vec3 operator+(float f_) const {
        return Vec3(x + f_, y + f_, z + f_);
    }

    /// lhs scalar cwise addition
    friend Vec3 operator+(const float d_, const Vec3 &vec) {
        return Vec3(d_ + vec.x, d_ + vec.y, d_ + vec.z);
    }

    Vec3 &operator+=(float f_) {
        x += f_;
        y += f_;
        z += f_;
        return *this;
    }

    // -------------
    // Substractions
    // -------------


    Vec3 operator-(const Vec3 &v_) const {
        return Vec3(x - v_.x, y - v_.y, z - v_.z);
    }


    Vec3 &operator-=(const Vec3 &v_) {
        x -= v_.x;
        y -= v_.y;
        z -= v_.z;
        return *this;
    }

    /// opposite vector

    Vec3 operator-() const {
        return Vec3(-x, -y, -z);
    }

    Vec3 operator-(float f_) const { return Vec3(x - f_, y - f_, z - f_); }

    /// lhs scalar cwise substraction
    friend Vec3 operator-(const float d_, const Vec3 &vec) {
        return Vec3(d_ - vec.x, d_ - vec.y, d_ - vec.z);
    }

    Vec3 &operator-=(float f_) {
        x -= f_;
        y -= f_;
        z -= f_;
        return *this;
    }

    // -------------
    // Comparisons
    // -------------

    bool operator!=(const Vec3 &v_) const {
        return (x != v_.x) | (y != v_.y) | (z != v_.z);
    }


    bool operator==(const Vec3 &d_) const {
        return (x == d_.x) && (y == d_.y) && (z == d_.z);
    }

    /// @note no mathematical meaning but useful to use Vec3 in std::map
    bool operator<(const Vec3 &v_) const {
        if (x != v_.x)
            return x < v_.x;
        else if (y != v_.y)
            return y < v_.y;
        else
            return z < v_.z;
    }

    // -------------
    // Divisions
    // -------------

    Vec3 operator/(const float d_) const {
        return Vec3(x / d_, y / d_, z / d_);
    }

    Vec3 &operator/=(const float d_) {
        x /= d_;
        y /= d_;
        z /= d_;
        return *this;
    }

    Vec3 operator/(const Vec3 &v_) const {
        return Vec3(x / v_.x, y / v_.y, z / v_.z);
    }

    // ----------------
    // Multiplication
    // ----------------

    /// rhs scalar multiplication
    Vec3 operator*(const float d_) const {
        return Vec3(x * d_, y * d_, z * d_);
    }

    /// lhs scalar multiplication
    friend Vec3 operator*(const float d_, const Vec3 &vec) {
        return Vec3(d_ * vec.x, d_ * vec.y, d_ * vec.z);
    }

    Vec3 &operator*=(const float d_) {
        x *= d_;
        y *= d_;
        z *= d_;
        return *this;
    }

    Vec3 operator*(const Vec3 &v_) const {
        return Vec3(x * v_.x, y * v_.y, z * v_.z);
    }

    Vec3 &operator*=(const Vec3 &d_) {
        x *= d_.x;
        y *= d_.y;
        z *= d_.z;
        return *this;
    }

    // -------------------------------------------------------------------------
    /// @name Operators on vector
    // -------------------------------------------------------------------------

    /// cross product
    Vec3 cross(const Vec3 &v_) const {
        return Vec3(y * v_.z - z * v_.y, z * v_.x - x * v_.z, x * v_.y - y * v_.x);
    }

    /// dot product
    float dot(const Vec3 &v_) const {
        return x * v_.x + y * v_.y + z * v_.z;
    }

    /// Compute the cotangente (i.e. 1./tan) between 'this' and v_
    float cotan(const Vec3 &v_) const {
        // cot(alpha ) = dot(v1, v2) / ||cross(v1, v2)||
        // = ||v1||*||v2||*cos( angle(v1, v2) ) / ||v1||*||v2|| * sin( angle(v1, v2) )
        // = cos( angle(v1, v2) ) / sin( angle(v1, v2) )
        // = 1 / tan( angle(v1, v2) )
        // = cot( angle(v1, v2) ) = cot( alpha )
        return (this->dot(v_)) / (this->cross(v_)).norm();
    }

    float norm_squared() const {
        return dot(*this);
    }

    Vec3 normalized() const {
        return (*this) * (float(1) / std::sqrt(norm_squared()));
    }

    float normalize() {
        float l = std::sqrt(norm_squared());
        float f = float(1) / l;
        x *= f;
        y *= f;
        z *= f;
        return l;
    }

    float safe_normalize(const float eps = 1e-10) {
        float l = std::sqrt(norm_squared());
        if (l > eps) {
            float f = float(1) / l;
            x *= f;
            y *= f;
            z *= f;
            return l;
        } else {
            x = float(1);
            y = float(0);
            z = float(0);
            return float(0);
        }
    }

    float norm() const {
        return std::sqrt(norm_squared());
    }

    // -------------------------------------------------------------------------
    /// @name Accessors
    // -------------------------------------------------------------------------

    const float &operator[](int i) const {
        return ((const float *) this)[i];
    }

    float &operator[](int i) {
        return ((float *) this)[i];
    }

    /// Conversion returns the memory address of the vector.
    /// Very convenient to pass a Vec pointer as a parameter to OpenGL:
    /// @code
    /// Vec3 pos, normal;
    /// glNormal3fv(normal);
    /// glVertex3fv(pos);
    /// @endcode
    operator const float *() const { return &x; }

    /// Conversion returns the memory address of the vector. (Non const version)
    operator float *() { return &x; }

};

enum Node_e {
    eUNIFORM,
    eOPEN_UNIFORM ///< Connected to the first and last control points
};

/**
 * @class Spline
 *
 * @brief Handling spline curves of arbitrary dimensions
 * @note This class use the efficient blossom algorithm to compute a position on
 * the curve.
 *
 * @tparam Point_t : type of a point operators such as '+' '*' must be correctly
 * overloaded. The default constructor must be defined to return the
 * null vector (0, 0 ,0 ...)
 * @tparam Real_t ; floating point reprensentation of the points
 * (float, double etc.)
 */
template<typename Point_t, typename Real_t>
class Spline {
public:

    /// Type of the nodal vector
    /// @param k : order of the spline (minimum is two)
    /// @param node_type : nodal vector type (uniform, open_uniform)
    /// This will define the behavior of the spline with its control points
    /// as well as its speed according to its parameter.
    Spline(int k = 2, Node_e node_type = eOPEN_UNIFORM);

    /// Set the position of the spline control points.
    void set_ctrl_points(const std::vector<Point_t> &point);

    /// Get the control points of the spline
    void get_ctrl_points(std::vector<Point_t> &points) const;

    /// The the nodal vector type
    void set_node_type(Node_e type);

    /// Evaluate position of the spline
    /// @param u : curve parameter ranging from [0; 1]
    Point_t eval_f(Real_t u) const;

    /// Evaluate speed of the spline
    Point_t eval_df(Real_t u) const;

    int get_order() const { return _k; }

private:
    // -------------------------------------------------------------------------
    /// @name Class tools
    // -------------------------------------------------------------------------

    void assert_splines() const;

    /// set value and size of the nodal vector depending on the current number
    /// of control points
    void set_nodal_vector();

    /// Set values of the nodal vector to be uniform
    void set_node_to_uniform();

    /// Set values of the nodal vector to be open uniform
    void set_node_to_open_uniform();

    /// Evaluate the equation of a splines using the blossom algorithm
    /// @param u : the curve parameter which range from the values
    /// [node[k-1]; node[point.size()]]
    /// @param point : the control points which size must be at least equal to
    /// the order of the spline (point.size() >= k)
    /// @param k : the spline order (degree == k-1)
    /// @param node : the nodal vector which defines the speed of the spline
    /// parameter u. The nodal vector size must be equal to (k + point.size())
    /// @param off : offset to apply to the nodal vector 'node' before reading
    /// from it. this is useful to compute derivatives.
    Point_t eval(Real_t u,
                 const std::vector<Point_t> &point,
                 int k,
                 const std::vector<Real_t> &node,
                 int off = 0) const;

    Point_t eval_rec(Real_t u,
                     std::vector<Point_t> p_in,
                     int k,
                     std::vector<Real_t> node_in) const;

    // -------------------------------------------------------------------------
    /// @name attributes
    // -------------------------------------------------------------------------

    Node_e _node_type;    ///< Nodal vector type
    int _k;                        ///< spline order
    std::vector<Point_t> _point;   ///< Control points
    std::vector<Point_t> _vec;     ///< Control points differences
    std::vector<Real_t> _node;    ///< Nodal vector
};


#include "splines.hpp"

#include <cassert>
#include <algorithm>

template<typename Point_t, typename Real_t>
Spline<Point_t, Real_t>::Spline(
        int k,
        Node_e node_type)
        :
        _node_type(node_type),
        _k(k),
        _point(_k),
        _vec(_k - 1),
        _node(_k + _point.size()) {
    assert_splines();
}

// -----------------------------------------------------------------------------

template<typename Point_t, typename Real_t>
void Spline<Point_t, Real_t>::set_ctrl_points(const std::vector<Point_t> &point) {
    _point = point;
    _vec.resize(_point.size() - 1);
    for (int i = 0; i < (int) _vec.size(); ++i)
        _vec[i] = _point[i + 1] - _point[i];
    set_nodal_vector();
    assert_splines();

    for (int i = 0; i < (int) _vec.size(); ++i)
        _vec[i] /= _node[_k + i] - _node[i + 1];
}

// -----------------------------------------------------------------------------

template<typename Point_t, typename Real_t>
void Spline<Point_t, Real_t>::get_ctrl_points(std::vector<Point_t> &points) const {
    points = _point;
}

// -----------------------------------------------------------------------------

/// The the nodal vector type
template<typename Point_t, typename Real_t>
void Spline<Point_t, Real_t>::set_node_type(Node_e type) {
    _node_type = type;
    set_nodal_vector();
    assert_splines();
}

// -----------------------------------------------------------------------------

template<typename Point_t, typename Real_t>
Point_t Spline<Point_t, Real_t>::eval_f(Real_t u) const {
    u = std::max(std::min(u, (Real_t) 1), (Real_t) 0); // clamp between [0 1]
    return eval(u, _point, _k, _node);
}

// -----------------------------------------------------------------------------

template<typename Point_t, typename Real_t>
Point_t Spline<Point_t, Real_t>::eval_df(Real_t u) const {
    u = std::max(std::min(u, (Real_t) 1), (Real_t) 0); // clamp between [0 1]
    return eval(u, _vec, (_k - 1), _node, 1) * (Real_t) (_k - 1);
}

// -----------------------------------------------------------------------------

template<typename Point_t, typename Real_t>
void Spline<Point_t, Real_t>::assert_splines() const {
    assert(_k > 1);
    assert((int) _point.size() >= _k);
    assert(_node.size() == (_k + _point.size()));
    assert(_point.size() == (_vec.size() + 1));
}

// -----------------------------------------------------------------------------

template<typename Point_t, typename Real_t>
void Spline<Point_t, Real_t>::set_nodal_vector() {
    if (_node_type == eOPEN_UNIFORM)
        set_node_to_open_uniform();
    else if (_node_type == eUNIFORM)
        set_node_to_uniform();
}

// -----------------------------------------------------------------------------

template<typename Point_t, typename Real_t>
void Spline<Point_t, Real_t>::set_node_to_uniform() {
    const int n = _point.size() - 1;
    _node.resize(_k + n + 1);

    Real_t step = (Real_t) 1 / (Real_t) (n - _k + 2);
    for (int i = 0; i < (int) _node.size(); ++i) {
        _node[i] = ((Real_t) i) * step - step * (Real_t) (_k - 1);
    }
}

// -----------------------------------------------------------------------------

template<typename Point_t, typename Real_t>
void Spline<Point_t, Real_t>::set_node_to_open_uniform() {
    _node.resize(_k + _point.size());

    int acc = 1;
    for (int i = 0; i < (int) _node.size(); ++i) {
        if (i < _k)
            _node[i] = 0.;
        else if (i >= ((int) _point.size() + 1))
            _node[i] = 1.;
        else {
            _node[i] = (Real_t) acc / (Real_t) (_point.size() + 1 - _k);
            acc++;
        }
    }
}

// -----------------------------------------------------------------------------

template<typename Point_t, typename Real_t>
Point_t Spline<Point_t, Real_t>::

eval(Real_t u,
     const std::vector<Point_t> &point,
     int k,
     const std::vector<Real_t> &node,
     int off) const {
    assert(k > 1);
    assert((int) point.size() >= k);
    assert_splines();
    int dec = 0;
    // TODO: better search with dychotomi ?
    // TODO: check for overflow
    while (u > node[dec + k + off])
        dec++;

    // TODO: use buffers in attributes for better performances ?
    std::vector<Point_t> p_rec(k, Point_t());
    for (int i = dec, j = 0; i < (dec + k); ++i, ++j)
        p_rec[j] = point[i];

    std::vector<Real_t> node_rec(k + k - 2, (Real_t) 0);
    for (int i = (dec + 1), j = 0; i < (dec + k + k - 1); ++i, ++j)
        node_rec[j] = node[i + off];

    return eval_rec(u, p_rec, k, node_rec);
}

// -----------------------------------------------------------------------------

template<typename Point_t, typename Real_t>
Point_t Spline<Point_t, Real_t>::

eval_rec(Real_t u,
         std::vector<Point_t> p_in,
         int k,
         std::vector<Real_t> node_in) const {
    if (p_in.size() == 1)
        return p_in[0];

    // TODO: use buffers in attributes for better performances ?
    std::vector<Point_t> p_out(k - 1, Point_t());
    for (int i = 0; i < (k - 1); ++i) {
        const Real_t n0 = node_in[i + k - 1];
        const Real_t n1 = node_in[i];
        const Real_t f0 = (n0 - u) / (n0 - n1);
        const Real_t f1 = (u - n1) / (n0 - n1);

        p_out[i] = p_in[i] * f0 + p_in[i + 1] * f1;
    }

    std::vector<Real_t> node_out(node_in.size() - 2);
    for (int i = 1, j = 0; i < ((int) node_in.size() - 1); ++i, ++j)
        node_out[j] = node_in[i];

    return eval_rec(u, p_out, (k - 1), node_out);
}


#endif // SPLINES_HPP
