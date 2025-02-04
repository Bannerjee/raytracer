#pragma once
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <random>

static constexpr double infinity = std::numeric_limits<double>::infinity();

double degrees_to_radians(double degrees) {
    return degrees * std::numbers::pi / 180.0;
}

double random_double(double min = 0,double max = 1)
{
    static std::mt19937_64 mt(std::random_device{}());
    std::uniform_real_distribution<double> dist(min,max);
    return dist(mt); 
}
struct interval 
{
    double min, max;

    interval() : min(+infinity), max(-infinity) {} 

    interval(double min, double max) : min(min), max(max) {}

    double size() const {
        return max - min;
    }

    bool contains(double x) const {
        return min <= x && x <= max;
    }

    bool surrounds(double x) const {
        return min < x && x < max;
    }

    double clamp(double x) const {
        if (x < min) return min;
        if (x > max) return max;
        return x;
    }

    static const interval empty, universe;
};

const interval interval::empty    = interval(+infinity, -infinity);
const interval interval::universe = interval(-infinity, +infinity);


