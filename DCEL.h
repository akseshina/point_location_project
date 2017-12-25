#pragma once

#include "point.h"
#include "double_arithmetics.h"
#include "geom_primitives.h"
#include "search_structure.h"

#include <vector>
#include <cstdio>
#include <stack>
#include <algorithm>
#include <map>
#include <cassert>
#include <unordered_map>


const int DEG_BOUND = 12;

typedef long double ld;


struct comp_segments {
    bool operator()(const Segment &lhs, const Segment &rhs) const {

        auto x1 = lhs.a->coord.x;
        auto y1 = lhs.a->coord.y;
        auto x2 = lhs.b->coord.x;
        auto y2 = lhs.b->coord.y;

        auto x3 = rhs.a->coord.x;
        auto y3 = rhs.a->coord.y;
        auto x4 = rhs.b->coord.x;
        auto y4 = rhs.b->coord.y;

        if (lhs.a == rhs.a && lhs.b == rhs.b)
            return false;

        auto yc = std::min(std::max(y1, y2), std::max(y3, y4));

        long double x_on_lhs = (yc - y1) * (long double) (x2 - x1) / (y2 - y1) + (long double) x1;
        long double x_on_rhs = x3;
        if (rhs.a != rhs.b)
            x_on_rhs = (yc - y3) * (long double) (x4 - x3) / (y4 - y3) + (long double) x3;

        return ls(x_on_lhs, x_on_rhs);
    }
};


void add_diagonal_to_list(Vertex *a, Vertex *b, std::vector<std::pair<Edge *, Edge *>> &diags) {
    if (a->v_id > b->v_id)
        std::swap(a, b);
    diags.push_back({a->one_starting_e->prev, b->one_starting_e});
    //std::cout << "new diagonal: " << a->v_id << ' ' << b->v_id << std::endl;
}


struct DCEL {
private:
    DCEL(const DCEL &);

    DCEL operator=(const DCEL &);

public:

    DCEL(const std::vector<Point> &poly_points){
        int N = poly_points.size();
        int max_v = N + 100, max_f = max_v * (DEG_BOUND + 5), max_e = max_f * 3;
        vertices = new Vertex[max_v]();
        faces = new Face[max_f]();
        edges = new Edge[max_e]();

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
            edges[i].adj_face = &(faces[1]);
            edges[i].e_id = i;

            // edge p[i + 1] -> p[i]
            int j = i + N;
            edges[j].prev = &(edges[(j + 1) % N + N]);
            edges[j].next = &(edges[(N + j - 1) % N + N]);
            edges[j].twin = &(edges[i]);
            edges[j].starting_v = &(vertices[(i + 1) % N]);
            edges[j].adj_face = &(faces[0]);
            edges[j].e_id = j;
        }

        faces[1].one_border_e = &(edges[0]);
        faces[1].f_id = 1;
        faces[0].one_border_e = &(edges[N]);
        faces[0].f_id = 0;
    }


    void build_triangulation_structure(int r, const std::vector<Point> &big_triangle) {
        int N = V;

        Edge *ein[3], *eout[3], *em_in[3], *em_out[3];

        for (int i = 0; i < 3; ++i) {
            ein[i] = &edges[E++];
            eout[i] = &edges[E++];
        }

        for (int i = 0; i < 3; ++i) {
            *ein[i] = Edge{ein[(i + 2) % 3], ein[(i + 1) % 3], eout[i], &vertices[N + i],
                           &faces[i == 0 ? 2 : 3], (int) (ein[i] - edges)};
            *eout[i] = Edge{eout[(i + 1) % 3], eout[(i + 2) % 3], ein[i], &vertices[N + (i + 1) % 3],
                            &faces[0], int (eout[i] - edges)};
        }

        for (int i = 0; i < 2; ++i) {
            em_in[i] = &edges[E++];
            em_out[i] = &edges[E++];
        }

        *em_in[0] = Edge{ein[2], &edges[2 * N - 1], em_out[0], &vertices[N], &faces[3],
                         (int) (em_in[0] - edges)};
        *em_out[0] = Edge{&edges[N], ein[0], em_in[0], &vertices[0], &faces[2],
                          (int) (em_out[0] - edges)};
        *em_in[1] = Edge{ein[0], &edges[N + (r + N - 1) % N], em_out[1], &vertices[N + 1], &faces[2],
                         (int) (em_in[1] - edges)};
        *em_out[1] = Edge{&edges[N + r], ein[1], em_in[1], &vertices[r], &faces[3],
                          (int) (em_out[1] - edges)};

        ein[2]->next = edges[2 * N - 1].prev = em_in[0];
        edges[N].next = ein[0]->prev = em_out[0];
        ein[0]->next = edges[N + (r + N - 1) % N].prev = em_in[1];
        edges[N + r].next = ein[1]->prev = em_out[1];

        for (int i = 0; i < 3; ++i)
            vertices[N + i] = Vertex{ein[i], big_triangle[i], N + i, v_type()};
        V += 3;

        F = 4;
        faces[0] = Face{eout[0], 0, 0};
        faces[1] = Face{&edges[0], 1, 1};
        faces[2] = Face{ein[0], 2, 0};
        faces[3] = Face{ein[1], 3, 0};

        for (int i = N; i < 2 * N; ++i)
            edges[i].adj_face = &faces[i < N + r ? 2 : 3];

        for (int i = 1; i <= 3; ++i)
            triangulate(&faces[i], i == 1);

        for (int i = 0; i < V; ++i)
            vertices[i].v_id = i;

        for (int i = 0; i < F; ++i) 
            assert((int) vertices_of_face(&faces[i]).size() == 3);
    }


    inline std::pair<int, int> deg(int v) {
        auto e0 = vertices[v].one_starting_e;
        if (!e0)
            return {0, 0};
        int cnt_all = 1, cnt = 0;
        auto e = e0->twin;
        for (; e->next != e0; e = e->next->twin, ++cnt_all)
            if (e->starting_v->v_id != -1)
                ++cnt;
        if (e->starting_v->v_id != -1)
            ++cnt;
        return {cnt_all, cnt};
    }


    void new_triangle(Vertex *a) { // insert diagonal b--c for triangle abc
        Edge *a_b = a->one_starting_e;
        Edge *c_a = a_b->prev;
        Vertex *b = a_b->next->starting_v;
        Vertex *c = c_a->starting_v;

        E += 2;
        edges[E - 2].e_id = E - 2;
        edges[E - 1].e_id = E - 1;
        Edge *b_c = &(edges[E - 2]);
        Edge *c_b = &(edges[E - 1]);

        Face *old_face = a_b->adj_face;
        old_face->one_border_e = c_b;

        F += 1;
        Face *new_face = &(faces[F - 1]);
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


    void add_diagonals(const std::vector<std::pair<Edge *, Edge *>> &diags, const std::vector<Vertex *> &vs) {
        // add new half-edges
        for (int i = 0; i < (int) diags.size(); ++i) {
            Edge *e_a_prev = diags[i].first;
            Vertex *a = e_a_prev->next->starting_v;

            Edge *e_b_next = diags[i].second;
            Vertex *b = e_b_next->starting_v;

            Vertex *a_next = e_a_prev->next->next->starting_v;
            // if a already has diag
            while (!(a_next->v_id > a->v_id && a_next->v_id < b->v_id)) {
                e_a_prev = e_a_prev->next->twin;
                a_next = e_a_prev->next->next->starting_v;
            }

            Vertex *b_prev = e_b_next->prev->starting_v;
            // if b already has diag
            while (!(b_prev->v_id > a->v_id && b_prev->v_id < b->v_id)) {
                e_b_next = e_b_next->prev->twin;
                b_prev = e_b_next->prev->starting_v;
            }

            Edge *e_a_next = e_a_prev->next;
            Edge *e_b_prev = e_b_next->prev;

            E += 2;
            edges[E - 2].e_id = E - 2;
            edges[E - 1].e_id = E - 1;
            Edge *ab = &(edges[E - 2]);
            Edge *ba = &(edges[E - 1]);

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
        Face *old_face = vs[0]->one_starting_e->adj_face;
        for (int j = 0; j < (int) vs.size(); ++j) {
            Edge *start_e = vs[j]->one_starting_e;
            if (start_e->adj_face == old_face) {

                F += 1;
                Face *new_face = &(faces[F - 1]);
                new_face->f_id = F - 1;
                new_face->one_border_e = start_e;
                start_e->adj_face = new_face;

                Edge *next_e = start_e->next;
                while (next_e != start_e) {
                    next_e->adj_face = new_face;
                    next_e = next_e->next;
                }
            }
        }

        // remove one useless face (swap last and old)
        Edge *start_e = faces[F - 1].one_border_e;
        start_e->adj_face = old_face;
        old_face->one_border_e = start_e;
        Edge *next_e = start_e->next;
        while (next_e != start_e) {
            next_e->adj_face = old_face;
            next_e = next_e->next;
        }
        --F;
    }


    void split_to_monotone(Face *f) {

        std::vector<Vertex *> vs = vertices_of_face(f);
        std::sort(vs.begin(), vs.end(), comp_y());
        assign_vertices_types(vs);

        int id_i = 0;
        vs[0]->v_id = id_i;
        Vertex *next_v = vs[0]->one_starting_e->next->starting_v;
        while (next_v != vs[0]) {
            ++id_i;
            next_v->v_id = id_i;
            next_v = next_v->one_starting_e->next->starting_v;
        }

        std::vector<std::pair<Edge *, Edge *>> diags;
        std::map<Segment, Vertex *, comp_segments> helper;

        for (int i = 0; i < (int) vs.size(); ++i) {
            Vertex *v_i = vs[i];
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
                        if (helper.find(e_i_prev) != helper.end() && helper[e_i_prev]->type == MERGE)
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


    void triangulate_monotone(Face *f) {

        std::vector<Vertex *> vs = vertices_of_face(f);
        std::sort(vs.begin(), vs.end(), comp_y());
        std::vector<Vertex *> s;
        s.push_back(vs[0]);
        s.push_back(vs[1]);

        for (int j = 2; j < int(vs.size()) - 1; ++j) {
            Vertex *e_end = s[0]->one_starting_e->next->starting_v;
            int orientation = 0;
            if (e_end == s[1]) {
                e_end = s[0]->one_starting_e->prev->starting_v;
                orientation = 1;
            }

            if (vs[j] == e_end) {
                for (int si = 0; si < (int) s.size() - 1; ++si)
                    new_triangle(s[si]);
                s.clear();
                s.push_back(vs[j - 1]);
                s.push_back(vs[j]);
            } else {
                Vertex *last = s.back();
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
            for (int si = 0; si < int(s.size()) - 2; ++si)
                new_triangle(s[si]);
        s.clear();
    }


    std::vector<Vertex *> vertices_of_face(Face *f) {
        // get vertices and reorient their one_starting_e
        Edge *start_e = f->one_border_e;
        std::vector<Vertex *> vs;
        Vertex *start_v = start_e->starting_v;
        vs.push_back(start_v);
        start_v->one_starting_e = f->one_border_e;
        Vertex *next_v = start_v->one_starting_e->next->starting_v;
        Vertex *prev_v = start_v;
        while (next_v != vs[0]) {
            next_v->one_starting_e = prev_v->one_starting_e->next;
            vs.push_back(next_v);
            prev_v = next_v;
            next_v = next_v->one_starting_e->next->starting_v;
        }
        return vs;
    }


    void assign_vertices_types(std::vector<Vertex *> & vs) {
        std::vector<v_type> res;
        for (int i = 0; i < (int) vs.size(); ++i) {
            Vertex *cur = vs[i];
            Vertex *prev = cur->one_starting_e->prev->starting_v;
            Vertex *next = cur->one_starting_e->next->starting_v;
            if (prev->coord.y < cur->coord.y && next->coord.y < cur->coord.y)
                if (left_rot(prev->coord, cur->coord, next->coord))
                    cur->type = START;
                else
                    cur->type = SPLIT;
            else if (prev->coord.y > cur->coord.y && next->coord.y > cur->coord.y)
                if (left_rot(prev->coord, cur->coord, next->coord))
                    cur->type = END;
                else
                    cur->type = MERGE;
            else
                cur->type = REGULAR;
        }
    }


    void kirkpatrick_build(int outerFace, Search_Structure &ss) {
        int real_v = V;
        std::vector<int> type(V); // outer = -2, deleted = -1, locked on this round = 1

        for (int i = 0; i < V; ++i)
            if (vertices[i].v_id == -1)
                type[i] = -1;
        {
            Edge *e = faces[outerFace].one_border_e;
            int u;
            while (!type[u = e->starting_v->v_id]) {
                assert(vertices[u].v_id == u);
                type[u] = -2, e = e->next;
            }
            assert(type[u] == -2);
        }

        while (real_v > 3) {
            int cnt = 0;
            for (auto &x: type)
                if (x == 1)
                    x = 0;

            std::pair<int, int> kp;
            for (int i = 0; i < V; ++i)
                if (!type[i] && (kp = deg(i)).second < DEG_BOUND) {
                    int k = kp.first;
                    auto &curv = vertices[i];
                    ++cnt;
                    type[i] = -1;
                    std::vector<Vertex *> polygon;
                    std::vector<bool> isInner;
                    std::vector<int> ids;
                    auto e0 = curv.one_starting_e, e = e0->twin;
                    for (int j = 0; j < k; ++j, e = e->next->twin) {
                        auto up = e->starting_v;
                        int u = up->v_id;
                        assert(u != -1);
                        if (!type[u])
                            type[u] = 1;
                        polygon.push_back(up);
                        isInner.push_back(e->adj_face->inner == 1);
                        ids.push_back(u);

                        up->one_starting_e = e->twin->next;
                        e->e_id = e->twin->e_id = -1;
                        e->adj_face->f_id = -1;
                    }

                    Face &new_f = faces[F];
                    new_f.one_border_e = polygon.front()->one_starting_e;
                    new_f.f_id = F;
                    ++F;

                    std::vector<Triangle_Part> old_tr, new_tr;

                    for (int j = 0; j < k; ++j) {
                        auto &p1 = polygon[j], &p2 = polygon[(j + 1) % k];
                        old_tr.emplace_back(std::vector<Vertex *>{&curv, p1, p2}, isInner[j]);
                        auto e1 = p1->one_starting_e, e2 = p2->one_starting_e;
                        e1->prev = e2, e2->next = e1;
                        e1->adj_face = &new_f;
                    }
                    curv.v_id = -1;

                    int F0 = F;
                    triangulate(&new_f, 0);
                    for (int j = 0; j < k; ++j)
                        polygon[j]->v_id = ids[j];

                    new_tr.emplace_back(vertices_of_face(&new_f), 0);
                    for (int j = F0; j < F; ++j)
                        new_tr.emplace_back(vertices_of_face(&faces[j]), 0);
                    for (auto tn: new_tr)
                        for (auto to: old_tr)
                            if (triangles_intersect(tn.triangle, to.triangle))
                                ss.add(tn, to), ss.root = ss.get_id(tn);
                }
            assert(cnt >= (real_v - 6) / 24);

            real_v -= cnt;
        }
    }


    void triangulate(Face *face, bool inner) {
        int F0 = F;
        split_to_monotone(face);
        int old_F = F;
        triangulate_monotone(face);
        for (int i = F0; i < old_F; ++i)
            triangulate_monotone(&faces[i]);
        if (inner) {
            face->inner = 1;
            for (int i = F0; i < F; ++i)
                faces[i].inner = 1;
        }
    }


    ~DCEL() {
        delete[] vertices;
        delete[] faces;
        delete[] edges;
    }

    Vertex *vertices;
    Face *faces;
    Edge *edges;
    int V, E, F;
};