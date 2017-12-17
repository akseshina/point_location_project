#pragma once

#include <vector>
#include <cstdio>
#include <stack>
#include <algorithm>
#include <map>

enum v_type {SPLIT, MERGE, START, END, REGULAR};

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
    v_type type;
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

struct Segment {
    Segment(Vertex * pa, Vertex * pb)
            : a(pa) , b(pb) {}
    Vertex * a;
    Vertex * b;
};

struct comp_y {
    bool operator() (const Vertex * a, const Vertex * b) {
        return a->coord.y > b->coord.y ||
               (a->coord.y == b->coord.y && a->coord.x < b->coord.x);
    }
};

struct comp_segments {
    bool operator() (const Segment & lhs, const Segment & rhs) const {
        int x1 = lhs.a->coord.x;
        int y1 = lhs.a->coord.y;
        int x2 = lhs.b->coord.x;
        int y2 = lhs.b->coord.y;
        int x0 = rhs.a->coord.x;
        int y0 = rhs.a->coord.y;
        double x_on_lhs = (y0 - y1) * (x2 - x1) / (y2 - y1) + x1;  // TODO: what if y0 == y1?
        return x_on_lhs < x0;
    }
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

        faces[1].one_border_e = &(edges[0]);
        faces[1].f_id = 1;
        faces[0].one_border_e = &(edges[N]);
        faces[0].f_id = 0;
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


    void add_diagonal_to_list(Vertex *a, Vertex *b, std::vector<std::pair<Edge *, Edge *>> &diags) {
        if (a->v_id > b->v_id)
            std::swap(a, b);
        diags.push_back({a->one_starting_e->prev, b->one_starting_e});
    }


    void add_diagonals(const std::vector<std::pair<Edge *, Edge *>> & diags, const std::vector<Vertex *> & vs) {

        size_t N = vs.size();

        // add new halph-edges
        for (int i = 0; i < diags.size(); ++i) {
            Edge * e_a_prev = diags[i].first;
            Vertex * a = e_a_prev->next->starting_v;
            // if a already has diag
            if (e_a_prev->next->next->starting_v->v_id != (a->v_id + 1) % N)
                e_a_prev = e_a_prev->next->twin;

            Edge * e_b_next = diags[i].second;
            Vertex * b = e_b_next->starting_v;
            // if b already has diag
            if (e_b_next->prev->starting_v->v_id != (N + b->v_id - 1) % N)
                e_b_next = e_b_next->prev->twin;

            Edge * e_a_next = e_a_prev->next;
            Edge * e_b_prev = e_b_next->prev;

            E += 2;
            edges[E - 2].e_id = E - 2;
            edges[E - 1].e_id = E - 1;
            Edge * ab = &(edges[E - 2]);
            Edge * ba = &(edges[E - 1]);

            ab->prev = e_a_prev;
            ab->next = e_b_next;
            ab->starting_v = a;
            ab->twin = ba;

            ba->prev = e_b_prev;
            ba->next = e_a_next;
            ba->starting_v = b;
            ba->twin = ab;

            e_a_prev->next = ab;
            e_a_next->prev = ba;
            e_b_prev->next = ba;
            e_b_next->prev = ab;
        }

        // add new faces
        Face * old_face = vs[0]->one_starting_e->adj_face;
        for (int j = 0; j < vs.size(); ++j) {
            Edge * start_e = vs[j]->one_starting_e;
            if (start_e->adj_face == old_face) {

                F += 1;
                Face * new_face = &(faces[F - 1]);
                new_face->f_id = F - 1;
                new_face->one_border_e = start_e;
                start_e->adj_face = new_face;

                Edge * next_e = start_e->next;
                while (next_e != start_e) {
                    next_e->adj_face = new_face;
                    next_e = next_e->next;
                }
            }
        }

        // remove one useless face (swap last and old)
        Edge * start_e = faces[F - 1].one_border_e;
        start_e->adj_face = old_face;
        old_face->one_border_e = start_e;
        Edge * next_e = start_e->next;
        while (next_e != start_e) {
            next_e->adj_face = old_face;
            next_e = next_e->next;
        }
        --F;
    }


    void split_to_monotone(Face * f) {

        std::vector<Vertex *> vs = vertices_of_face(f);
        std::sort(vs.begin(), vs.end(), comp_y());
        assign_vertices_types(vs);

        int id_i = 0;
        vs[0]->v_id = id_i;
        Vertex * next_v = vs[0]->one_starting_e->next->starting_v;
        while (next_v != vs[0]) {
            ++id_i;
            next_v->v_id = id_i;
            next_v = next_v->one_starting_e->next->starting_v;
        }
        
        std::vector<std::pair<Edge *, Edge *>> diags;
        std::map<Segment, Vertex *, comp_segments> helper;

        for (int i = 0; i < vs.size(); ++i) {
            Vertex * v_i = vs[i];
            Segment e_i(v_i, v_i->one_starting_e->next->starting_v);
            Segment e_i_prev(v_i->one_starting_e->prev->starting_v, v_i);

            switch (v_i->type) {
                case START: {
                    // Insert e_i in T, helper(e_i) <- v_i
                    helper[e_i] = v_i;
                    break;
                }

                case END: {
                    // if (Type_of_vertex(helper(e_i_prev) = 'merge')
                    //     Insert edge(v_i, helper(e_i_prev)) in D
                    // Delete e_i_prev from T
                    if (helper[e_i_prev]->type == MERGE)
                        add_diagonal_to_list(v_i, helper[e_i_prev], diags);
                    helper.erase(e_i_prev);
                    break;
                }

                case SPLIT: {
                    // edge e_j = l intersect T
                    // Search e_j in T
                    Segment e_j = std::prev(helper.lower_bound({v_i, v_i}))->first;

                    // Insert edge(v_i, helper(e_j)) in D
                    add_diagonal_to_list(v_i, helper[e_j], diags);

                    // helper(e_j) <- v_i
                    helper[e_j] = v_i;

                    // Insert e_i in T, helper(e_i) <- v_i
                    helper[e_i] = v_i;
                    break;
                }

                case MERGE: {
                    // if (Type_of_vertex(helper(e_i_prev) = 'merge')
                    //     Insert edge(v_i, helper(e_i_prev)) in D
                    // Delete e_i_prev from T
                    if (helper[e_i_prev]->type == MERGE)
                        add_diagonal_to_list(v_i, helper[e_i_prev], diags);
                    helper.erase(e_i_prev);

                    // edge e_j = l intersect P
                    // Search e_j in T
                    Segment e_j = std::prev(helper.lower_bound({v_i, v_i}))->first;

                    // if (Type_of_vertex(helper(e_j) = 'merge')
                    //     Insert edge(v_i, helper(e_j)) in D
                    if (helper[e_j]->type == MERGE)
                        add_diagonal_to_list(v_i, helper[e_j], diags);

                    // helper(e_j) <- v_i
                    helper[e_j] = v_i;
                    break;
                }

                case REGULAR: {
                    // if (interior of P lies to the right of v_i)
                    if (e_i.b->coord.y < e_i_prev.a->coord.y) {
                        // if (Type_of_vertex(helper(e_i_prev) = 'merge')
                        //     Insert edge(v_i, helper(e_i_prev)) in D
                        // Delete e_i_prev from T
                        if (helper[e_i_prev]->type == MERGE)
                            add_diagonal_to_list(v_i, helper[e_i_prev], diags);
                        helper.erase(e_i_prev);
                        // Insert e_i in T, helper(e_i) <- v_i
                        helper[e_i] = v_i;

                    } else {
                        // edge e_j = l intersect P
                        // Search e_j in T
                        Segment e_j = std::prev(helper.lower_bound({v_i, v_i}))->first;

                        // if (Type_of_vertex(helper(e_j) = 'merge')
                        //     Insert edge(v_i, helper(e_j)) in D
                        if (helper[e_j]->type == MERGE)
                            add_diagonal_to_list(v_i, helper[e_j], diags);

                        // helper(e_j) <- v_i
                        helper[e_j] = v_i;
                    }
                    break;
                }
            }
        }

        add_diagonals(diags, vs);
    }


    void triangulate_monotone(Face * f) {

        std::vector<Vertex *> vs = vertices_of_face(f);
        std::sort(vs.begin(), vs.end(), comp_y());
        std::vector<Vertex *> s;
        s.push_back(vs[0]);
        s.push_back(vs[1]);

        for (int j = 2; j < vs.size() - 1; ++j) {
            Vertex * e_end = s[0]->one_starting_e->next->starting_v;
            int orientation = 0;
            if (e_end == s[1]) {
                e_end = s[0]->one_starting_e->prev->starting_v;
                orientation = 1;
            }

            if (vs[j] == e_end) { // Текущая вершина v_j является нижним концом стороны e,
                                  // ограничивающей воронку
                for (int si = 0; si < s.size() - 1; ++si)
                    new_triangle(s[si]);
                s.clear();
                s.push_back(vs[j - 1]);
                s.push_back(vs[j]);
            }

            else {                // Вершина v_j принадлежит последовательной цепи вершин, добавленных в S
                Vertex * last = s.back();
                s.pop_back();
                if (orientation == 0)
                    while (!s.empty() && left_rot(vs[j]->coord, last->coord, s.back()->coord)) {
                        new_triangle(last);
                        last = s.back();
                        s.pop_back();
                    }
                else
                    while (!s.empty() && right_rot(vs[j]->coord, last->coord, s.back()->coord)) {
                        new_triangle(last);
                        last = s.back();
                        s.pop_back();
                    }
                s.push_back(last);
                s.push_back(vs[j]);
            }
        }

        if (s.size() > 2)
            for (int si = 0; si < s.size() - 1; ++si)
                new_triangle(s[si]);
        s.clear();
    }


    std::vector<Vertex*> vertices_of_face(Face * f) {
        // get vertices and reorient their one_starting_e
        Edge * start_e = f->one_border_e;
        std::vector<Vertex *> vs;
        Vertex * start_v = start_e->starting_v;
        vs.push_back(start_v);
        start_v->one_starting_e = f->one_border_e;
        Vertex * next_v = start_v->one_starting_e->next->starting_v;
        Vertex * prev_v = start_v;
        while (next_v != vs[0]) {
            next_v->one_starting_e = prev_v->one_starting_e->next;
            vs.push_back(next_v);
            prev_v = next_v;
            next_v = next_v->one_starting_e->next->starting_v;
        }
        return vs;
    }


    void assign_vertices_types(std::vector<Vertex *> vs) {
        std::vector<v_type> res;
        for (int i = 0; i < vs.size(); ++i) {
            Vertex * cur = vs[i];
            Vertex * prev = cur->one_starting_e->prev->starting_v;
            Vertex * next = cur->one_starting_e->next->starting_v;
            if (prev->coord.y < cur->coord.y && next->coord.y < cur->coord.y)
                if (left_rot(prev->coord, cur->coord, next->coord))
                    cur->type = START;
                else
                    cur->type = SPLIT;
            else
                if (prev->coord.y > cur->coord.y && next->coord.y > cur->coord.y)
                    if (left_rot(prev->coord, cur->coord, next->coord))
                        cur->type = END;
                    else
                        cur->type = MERGE;
                else
                    cur->type = REGULAR;
        }
    }


    Vertex vertices[MAX_V];
    Face faces[MAX_F];
    Edge edges[MAX_E];
    int V, E, F;

};



