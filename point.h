#pragma once

#include "double_arithmetics.h"

typedef long double ld;

struct Point {
    Point()
            : x(0), y(0) {}

    Point(ld cur_x, ld cur_y)
            : x(cur_x), y(cur_y) {}

    ld x, y;

    inline Point operator-(const Point &b) const {
        return Point(x - b.x, y - b.y);
    }

    inline ld operator*(const Point &b) const {
        return x * b.x + y * b.y;
    }

    inline ld operator%(const Point &b) const {
        return x * b.y - y * b.x;
    }

    inline bool operator<(const Point &b) const {
        return ls(x, b.x) || (eq(x, b.x) && ls(y, b.y));
    }

    inline bool operator==(const Point &b) const {
        return eq(x, b.x) && eq(y, b.y);
    }

    inline void rotate(ld s, ld c){
        ld x0 = c * x - s * y, y0 = c * y + s * x;
        x = x0, y = y0;
    }
};

