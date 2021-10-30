#ifndef __EXTRAFUNCS_INCLUDED__
#define __EXTRAFUNCS_INCLUDED__

#include <math.h>
#include <cmath>
#include <sstream>
#include <string>    // string
#include <vector>    // vector
#include <iostream>

inline int pmod(int m, int n) {
    return (m % n + n) % n;
}

inline double fpmod(double m, double n) {
    return fmod((fmod(m,n)+n),n);
}

inline int sgn(double x) {
    return (x > 0) - (x < 0);
}

inline int ceildiv(int x, int y) {
    return (x + y - 1) / y;
}

inline int smaller(int x, int y) {
    if (x<y) {
        return x;
    }
    return y;
}

inline double smaller(double x, double y) {
    if (x<y) {
        return x;
    }
    return y;
}

inline int larger(int x, int y) {
    if (x>y) {
        return x;
    }
    return y;
}

inline double larger(double x, double y) {
    if (x>y) {
        return x;
    }
    return y;
}

inline int int_log10( int val )
{
    int num = 0;
    val/=10;
    while (val>0) {
        num++;
        val/=10;
    }
    return num;
}

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

template <typename T>
bool contained_in_Vec(const std::vector<T> & vec,T val) {
    for (unsigned i=0;i<vec.size();i++) {
        if (vec[i] == val) {
            return true;
        }
    }
    return false;
}

template <typename T>
void print_vector(const std::vector<T> & vec) {
    for (unsigned i=0;i<vec.size();i++) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}


std::vector<double> strVec2doubleVec(const std::vector<std::string> & strvec);
std::vector<int> strVec2intVec(const std::vector<std::string> & strvec);





#endif
