#pragma once

#include "point.h"

#include <cassert>

const int MAX_V = 2e6 + 10;
const int MAX_E = MAX_V;
const int MAX_F = MAX_V;

inline int sign(ld a) {
    return (a > 0) - (a < 0);
}

enum v_type {
    SPLIT, MERGE, START, END, REGULAR
};


struct Edge;


struct Vertex {
    Edge *one_starting_e;  /* any half-edge which starts at this Vertex */
    Point coord;
    int v_id;
    v_type type;
};


struct Face {
    Edge *one_border_e;  /* any half-edge on the border */
    int f_id;
    int inner;
};


struct Edge {
    Edge *prev;
    Edge *next;
    Edge *twin;
    Vertex *starting_v;
    Face *adj_face;
    int e_id;
};


struct Segment {
    Segment(Vertex *pa, Vertex *pb)
            : a(pa), b(pb) {}

    Vertex *a;
    Vertex *b;
};


struct comp_y {
    bool operator()(const Vertex *a, const Vertex *b) {
        return gr(a->coord.y, b->coord.y) ||
               (eq(a->coord.y, b->coord.y) && ls(a->coord.x, b->coord.x));
    }
};


struct Triangle {
    Point v[3];

    Triangle() {}
    Triangle(const std::vector<Point> &a) {
        assert((int) a.size() == 3);
        for (int k = 0; k < 3; ++k)
            v[k] = a[k];
        if ((v[1] - v[0]) % (v[2] - v[0]) < 0)
            std::swap(v[1], v[2]);
    }

    bool contains(const Point &p) {
        return geq((v[1] - v[0]) % (p - v[0]), 0) &&
               leq((v[2] - v[0]) % (p - v[0]), 0) &&
               geq((v[2] - v[1]) % (p - v[1]), 0);
    }
};

struct Triangle_Part {
    std::vector<int> nums;
    Triangle triangle;
    bool inner;
    Triangle_Part(std::vector<Vertex*> v, bool in):inner(in){
        assert((int)v.size() == 3);
        triangle = Triangle({v[0]->coord, v[1]->coord, v[2]->coord});
        nums = {v[0]->v_id, v[1]->v_id, v[2]->v_id};
        std::sort(nums.begin(), nums.end());
    }
};


inline bool on_interval(const Point &p, const Point &a, const Point &b) {
    return eq((p - a) % (b - a), 0) && geq((p - a) * (b - a), 0) && geq((p - b) * (a - b), 0);
}

bool intervals_intersect(const Point &a, const Point &b, const Point &c, const Point &d) {
    if (on_interval(a, c, d) || on_interval(b, c, d) || on_interval(c, a, b) || on_interval(d, a, b))
        return 1;
    return sign((a - b) % (c - b)) * sign((a - b) % (d - b)) < 0 &&
           sign((c - d) % (a - d)) * sign((c - d) % (b - d)) < 0;
}

bool triangles_intersect(const Triangle &a, const Triangle &b) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (intervals_intersect(a.v[i], a.v[(i + 1) % 3], b.v[j], b.v[(j + 1) % 3]))
                return 1;
    return 0;
}


ld rot_matrix(Point a, Point b, Point c) {
    return a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y);
}

bool left_rot(Point a, Point b, Point c) {
    return gr(rot_matrix(a, b, c), 0);
}

bool right_rot(Point a, Point b, Point c) {
    return ls(rot_matrix(a, b, c), 0);
}