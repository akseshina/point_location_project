#pragma once

#include <algorithm>

typedef long double ld;

const ld EPS = 1e-9;


inline bool gr(ld a, ld b){
    return a > b + EPS;
}

inline bool ls(ld a, ld b){
    return a < b - EPS;
}

inline bool eq(ld a, ld b){
    return std::abs(a - b) < EPS;
}

inline bool geq(ld a, ld b){
    return a >= b - EPS;
}

inline bool leq(ld a, ld b){
    return a <= b + EPS;
}
