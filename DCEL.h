#pragma once

#include <vector>
#include <cstdio>
#include <stack>
#include <algorithm>

//TODO: calculate real parameters
const int MAX_V = 1000;
const int MAX_E = 1000;
const int MAX_F = 1000;

struct Edge;

struct point {
    point()
            : x(0), y(0) {}

    point(int cur_x, int cur_y)
            : x(cur_x), y(cur_y) {}

    int x, y;
};

double rot_matrix(point a, point b, point c) {
    return a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y);
}

bool left_rot(point a, point b, point c) {
    return rot_matrix(a, b, c) > 0;
}

bool right_rot(point a, point b, point c) {
    return rot_matrix(a, b, c) < 0;
}

struct Vertex {
    Edge * one_starting_e;  /* any half-edge which starts at this Vertex */
    point coord;
    int v_id;
};

struct Face {
    Edge * one_border_e;  /* any half-edge on the border */
    int f_id;
};

struct Edge {
    Edge * prev;
    Edge * next;
    Edge * twin;
    Vertex * starting_v;
    Face * adj_face;
    int e_id;
};


struct DCEL {

    DCEL(const std::vector<point> & poly_points) {
        //TODO: reverse points if they are not CCW

        unsigned long N = poly_points.size();
        V = N;
        E = 2 * N;  // first N are for inner cycle, second N are for outer
        F = 2;      // 1 is for inner, 0 is for outer

        for (int i = 0; i < N; ++i) {
            vertices[i].coord = poly_points[i];
            vertices[i].one_starting_e = &(edges[i]);
            vertices[i].v_id = i;

            // edge p[i] -> p[i + 1]
            edges[i].prev = &(edges[(N + i - 1) % N]);
            edges[i].next = &(edges[(i + 1) % N]);
            edges[i].twin = &(edges[i + N]);
            edges[i].starting_v = &(vertices[i]);
            edges[i].adj_face = &(faces[0]);
            edges[i].e_id = i;

            // edge p[i + 1] -> p[i]
            int j = i + N;
            edges[j].prev = &(edges[(j + 1) % N + N]);
            edges[j].next = &(edges[(N + j - 1) % N + N]);
            edges[j].twin = &(edges[i]);
            edges[j].starting_v = &(vertices[(i + 1) % N]);
            edges[j].adj_face = &(faces[1]);
            edges[j].e_id = j;
        }

        faces[0].one_border_e = &(edges[0]);
        faces[0].f_id = 1;
        faces[1].one_border_e = &(edges[N]);
        faces[1].f_id = 0;
    }

    void new_triangle(Vertex * a) { // insert diagonal b--c for triangle abc

        Edge * a_b = a->one_starting_e;
        Edge * c_a = a_b->prev;
        Vertex * b = a_b->next->starting_v;
        Vertex * c = c_a->starting_v;

        E += 2;
        edges[E - 2].e_id = E - 2;
        edges[E - 1].e_id = E - 1;
        Edge * b_c = &(edges[E - 2]);
        Edge * c_b = &(edges[E - 1]);

        Face * old_face = a_b->adj_face;
        old_face->one_border_e = c_b;

        F += 1;
        Face * new_face = &(faces[F - 1]);
        new_face->f_id = F - 1;
        new_face->one_border_e = b_c;
        a_b->adj_face = new_face;
        c_a->adj_face = new_face;

        b_c->prev = a_b;
        b_c->next = c_a;
        b_c->adj_face = new_face;
        b_c->twin = c_b;
        b_c->starting_v = b;

        c_b->prev = c_a->prev;
        c_b->next = a_b->next;
        c_b->adj_face = old_face;
        c_b->twin = b_c;
        c_b->starting_v = c;

        a_b->next->prev = c_b;
        a_b->next = b_c;

        c_a->prev->next = c_b;
        c_a->prev = b_c;

        c->one_starting_e = c_b;
    }

    struct less_by_y_key {
        inline bool operator() (const Vertex * a, const Vertex * b) {
            return a->coord.y > b->coord.y ||
                    (a->coord.y == b->coord.y && a->coord.x < b->coord.x);
        }
    };

    void TriangulateMonotonePolygon(std::vector<Vertex *> vs) {
        // TODO: get face, reorient one_starting_e for vertices

        std::sort(vs.begin(), vs.end(), less_by_y_key());
        std::vector<Vertex *> s;
        s.push_back(vs[0]);
        s.push_back(vs[1]);
        for (int j = 2; j < vs.size() - 1; ++j) {
            Vertex * e_end = s[0]->one_starting_e->next->starting_v;

            if (vs[j] == e_end) { // Текущая вершина v_j является нижним концом стороны e,
                                  // ограничивающего воронку
                for (int si = 0; si < s.size() - 1; ++si)
                    new_triangle(s[si]);
                s.clear();
                s.push_back(vs[j]);
                s.push_back(vs[j - 1]);
            }

            else {                // Вершина v_j принадлежит последовательной цепи вершин, добавленных в S
                Vertex * last = s.back();
                s.pop_back();
                while (left_rot(vs[j]->coord, last->coord, s.back()->coord)) {
                    new_triangle(last);
                    last = s.back();
                    s.pop_back();
                }
                s.push_back(last);
                s.push_back(vs[j]);
            }
        }

        if (vs.size() > 2)
            for (int si = 0; si < s.size() - 1; ++si)
                new_triangle(s[si]);
        s.clear();
    }


    Vertex vertices[MAX_V];
    Face faces[MAX_F];
    Edge edges[MAX_E];
    int V, E, F;
};



