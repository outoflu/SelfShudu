#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <iostream>
#include <cmath>
#include <cassert>
#include <initializer_list>

template <int n>
struct vec
{
    double data[n];
    vec(std::initializer_list<double> &&_list)
    {
        static_assert(std::is_integral<int>::value && n > 0);
        auto list = std::forward<std::initializer_list<double>>(_list);
        int size = list.size();
        auto iter = list.begin();
        for (int i = 0; i < n; i++, iter++)
        {
            data[i] = i > size ? 0 : *iter;
        }
    }
    vec()
    {
        for (int i = 0; i < n; i++)
        {
            data[i] = 0;
        }
    }
    double &operator[](const int x)
    {
        assert(x >= 0 && x < n);
        return data[x];
    }
    double operator[](const int x) const
    {
        assert(x >= 0 && x < n);
        return data[x];
    }
    double norm2()
    {
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
            sum += i * i;
        }
        return sum;
    };
    double norm()
    {
        return sqrt(norm2());
    };
};
// 数乘
template <int n>
vec<n> operator*(vec<n> &lhs, double rhs)
{
    vec<n> ans;
    for (int i = 0; i < n; i++)
    {
        ans[i] = lhs[i] * rhs;
    }
    return std::move(ans);
}

template <int n>
vec<n> operator*(double lhs, vec<n> &rhs)
{
    vec<n> res;
    for (int i = 0; i < n; i++)
    {
        res[i] = lhs * rhs[i];
    }
    return std::move(res);
}

template <int n>
vec<n> operator/(vec<n> &lhs, double rhs)
{
    vec<n> res;
    for (int i = 0; i < n; i++)
    {
        res[i] = lhs[i] / rhs;
    }
    return std::move(res);
}
/// 点乘
template <int n>
double operator*(vec<n> &lhs, vec<n> &rhs)
{
    double res = 0;
    for (int i = 0; i < n; i++)
    {
        res += lhs[i] * rhs[i];
    }
    return res;
}

template <int n>
vec<n> operator-(vec<n> &lhs, vec<n> &rhs)
{
    vec<n> res;
    for (int i = 0; i < n; i++)
    {
        res[i] = lhs[i] - rhs[i];
    }
    return res;
}

template <int n>
vec<n> operator+(vec<n> &lhs, vec<n> &rhs)
{
    vec<n> res;
    for (int i = 0; i < n; i++)
    {
        res[i] = lhs[i] + rhs[i];
    }
    return res;
}

template <int n>
std::ostream &operator<<(std::ostream &os, vec<n> &v)
{
    for (int i = 0; i < n; i++)
    {
        os << " " << v[i];
    }
    os << std::endl;
    return os;
}

/// specify vec<2>
template <>
struct vec<2>
{
    double x, y;
    vec(double _x, double _y) : x(_x), y(_y){};
    
    double norm2()
    {
        return x * x + y * y;
    }
    double norm()
    {
        return std::sqrt(norm2());
    }

    double operator[](int idx) const{
        assert(idx >= 0 && idx < 2);
        return idx == 0? x : y;
    }
};

template <>
struct vec<3>
{
    double x, y, z;
    vec(double _x, double _y, double _z) : x(_x), y(_y), z(_z)
    {
    }
    vec()
    {
    }
    double norm2()
    {
        return x * x + y * y + z * z;
    }
    double norm()
    {
        return std::sqrt(norm2());
    }
    double &operator[] (int idx){
        assert(idx >= 0 && idx < 3);
        return idx == 0? x : idx == 1? y : z;
    }
    double operator[](int idx) const{
        assert(idx >= 0 && idx < 3);
        return idx == 0? x : idx == 1? y : z;
    }
};

vec<3> cross(vec<3> &lhs, vec<3> &rhs)
{   
    return vec<3>{lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x};
}


template <int rows, int cols>
struct matrix
{
    vec<cols> data[rows];
    vec<cols> &operator[](int idx) { return data[idx]; }
    const vec<cols> &operator[](int idx) const { return data[idx]; }

    vec<rows> col(const int idx) const
    {
        vec<rows> ret;
        for (int i = rows; i--; ret[i] = data[i][idx])
            ;
        return ret;
    }

    void set_col(const int idx, const vec<rows> &v)
    {
        for (int i = rows; i--; data[i][idx] = v[i])
            ;
    }

    static matrix<rows, cols> unit_matrix()
    {
        matrix<rows, cols> ret;
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                ret[i][j] = (i == j);
            }
        }
        return ret;
    }
    // 余子式
    matrix<rows - 1, cols - 1> get_minor(int row_index, int col_index)
    {
        matrix<rows - 1, cols - 1> ret;
        for (int i = 0; i < rows - 1; i++)
        {
            for (int j = 0; j < cols - 1; j++)
            {
                ret[i][j] = data[i < row_index ? i : i + 1][j < col_index ? j : j + 1];
            }
        }
        return ret;
    }
    double det()
    {
        return dt<rows>::det(*this);
    }
    double cofactor(const int row_index, const int col_index)
    {
        return get_minor(row_index, col_index).det() * ((i + j) % 2 ? -1 : 1);
    }
    /// @brief 余子式矩阵
    /// @return 
    matrix<rows, cols> minor_matrix()
    {
        matrix<rows, cols> ret;
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                ret[i][j] = cofactor(i, j);
            }
        }
        return ret;
    }
    /// @brief 余子式阵的转置矩阵才是伴随矩阵
    /// @return 
    matrix<rows, cols> companion_matrix(){
        this->minor_matrix().transpose_matrix();
    }
    /// @brief 求逆矩阵第一步
    /// @return 
    matrix<rows,cols> inverse_trans(){
        matrix<rows,cols> ret = minor_matrix();
        double _div=ret[0][0]*data[0][0];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                ret[i][j] /= _div;
            }
        }
        return ret;
    }
    /// @brief 逆矩阵
    /// @return 
    matrix<rows,cols> inverse_matrix(){
        return this->inverse_trans().transpose_matrix();
    }
    /// @brief 转置矩阵
    /// @return 
    matrix<rows, cols> transpose_matrix(){

        matrix<cols,rows> ret;
        for (int i = 0; i < rows; i++){
            ret[i]=this->col(i);
        }
        return ret;
    }
};
template<int rows,int cols>
matrix<rows, cols> operator*(const matrix<rows, cols>& lhs,double rhs) {
    matrix<rows, cols> ret;
        for (int i = 0; i < rows; i++) {
            ret[i]=lhs[i]*rhs;
        }
        return ret;
} 

template<int rows,int cols>
vec<rows> operator*(const matrix<rows, cols>& lhs,vec<cols>& rhs) {
    vec<rows> ret;
    for (int i = 0; i < rows; i++){
        ret[i]=lhs[i]*rhs;
    }
    return ret;
}

template<int r1,int r2,int c2>
matrix<r1,c2> operator*(const matrix<r1,r2> &lhs,const matrix<r2,c2> &rhs){
    matrix<r1,c2> ret;
    for (int i = 0; i < r1; i++){
        for (int j = 0; j < c2; j++){
            ret[i][j]=lhs[i]*rhs.col(j);
        }
    }
    return ret;
}

template<int rows,int cols>
matrix<rows,cols> operator/ (const matrix<rows,cols> &lhs,double rhs){
    matrix<rows,cols> ret;
    for (int i=0;i<rows;i++){
        ret[i]=lhs[i]/rhs;
    }
    return ret;
}

template <int rows,int cols>
matrix<rows,cols> operator+ (const matrix<rows,cols> &lhs,const matrix<rows,cols> &rhs){
    matrix<rows,cols> ret;
    for (int i=0;i<rows;i++){
        ret[i]=lhs[i]+rhs[i];
    }
    return ret;
}

template <int rows,int cols>
matrix<rows,cols> operator- (const matrix<rows,cols> &lhs,const matrix<rows,cols> &rhs){
    matrix<rows,cols> ret;
    for (int i=0;i<rows;i++){
        ret[i]=lhs[i]-rhs[i];
    }
    return ret;
}

template <int rows,int cols>
std::ostream &operator<<(std::ostream &os, const matrix<rows, cols> &m)
{
    for (int i=0;i<rows;++i) os<<m[i]; 
}

template <int rows>
struct dt
{
    static double det(const matrix<rows, cols> &src)
    {
        double ret = 0;
        for (int i = 0; i < rows; i++)
        {
            ret += src[0][i] * src.cofactor(0, i);
        }
        return ret;
    }
};

template<>
struct dt<3>{
    //由于重载了vec<2> vec<3> 这里应该写上
    static double det(const matrix<3, 3> &src){
        double ret = 0;
        ret+= src[0].x*src[1].y*src[2].z+src[0].y*src[1].z*src[2].x+src[0].z*src[1].x*src[2].y;
        ret-= src[2].x*src[1].y*src[0].z-src[2].y*src[1].z*src[0].x-src[2].z*src[1].x*src[0].y;
        return ret;
    }
};
template<>
struct dt<2>{
    //由于重载了vec<2> vec<3> 这里应该写上
    static double det(const matrix<2,2> &src){
        return src[0].x*src[1].y-src[0].y*src[1].x;
    }
};
template <>
struct dt<1>
{
    static double det(const matrix<1, 1> &src)
    {
        return src[0][0];
    }
};
#endif
