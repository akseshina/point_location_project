#include <iostream>
#include "DCEL.h"
#include "read_file.h"
#include <unordered_map>

/*
В первой строке дано целое число T — количество тестов.
Далее идут T тестов. Тесты разделены переводом строки.
N
Далее N точек — вершины многоугольника.
K
Далее K точек — запросы.
Все координаты — целые числа по модулю не превосходящие 10^9
*/


void make_ccw(std::vector<point> &points) {
    std::rotate(points.begin(), std::min_element(points.begin(), points.end()), points.end());
    const point &a = points[1], &b = points[0], &c = points.back();
    if ((c - b) % (a - b) > 0)
        reverse(points.begin() + 1, points.end());
}

const long long INF = 1e9 + 10;
const point LEFT(-2 * INF, -INF - 1), RIGHT(2 * INF, -INF), TOP(0, 3 * INF);
const point bigTriangle[3] = {LEFT, RIGHT, TOP};

void triangulate(DCEL &dcel, Face *face, bool inner) {
	int F0 = dcel.F;
    dcel.split_to_monotone(face);
    int old_F = dcel.F;
    dcel.triangulate_monotone(face);
    for (int i = F0; i < old_F; ++i)
        dcel.triangulate_monotone(&dcel.faces[i]);
    if(inner){
	    face->inner = 1;
        for (int i = F0; i < dcel.F; ++i)
            dcel.faces[i].inner = 1;
    }
}

std::pair<ld, ld> rt(std::vector<point> &points) {
    ld s = rand() % 100 / 100.0, c = sqrtl(1 - s * s);
//    ld s = 0, c = 1;
    for (auto &p: points)
        p.rt(s, c);
    return {s, c};
}

int main() {
    //freopen("input.txt", "r", stdin);

    int T, N, K;
    scanf("%d", &T);

    for (int t = 0; t < T; ++t) {
        scanf("%d", &N);
        assert(N >= 3);
        std::vector<point> poly_points = read_points(N);
        auto sc = rt(poly_points);
/*        for(auto &p: poly_points)
            printf("%.3f %.3f\n", (double)p.x, (double)p.y);*/

        make_ccw(poly_points); // leftest bottoms point now has index 0
        DCEL dcel(poly_points);
		int r = max_element(poly_points.begin(), poly_points.end()) - poly_points.begin();
		assert(r > 0);
		Edge *ein[3], *eout[3], *em_in[3], *em_out[3];

        for(int i = 0; i < 3; ++i){
        	ein[i] = &dcel.edges[dcel.E++];
        	eout[i] = &dcel.edges[dcel.E++];
        }
        for(int i = 0; i < 3; ++i){
		    *ein[i] = Edge{ein[(i + 2) % 3], ein[(i + 1) % 3], eout[i], &dcel.vertices[N + i], &dcel.faces[i == 0 ? 2 : 3], ein[i] - dcel.edges};
		    *eout[i] = Edge{eout[(i + 1) % 3], eout[(i + 2) % 3], ein[i], &dcel.vertices[N + (i + 1) % 3], &dcel.faces[0], eout[i] - dcel.edges};
        }
        for(int i = 0; i < 2; ++i){
	        em_in[i] = &dcel.edges[dcel.E++];
	    	em_out[i] = &dcel.edges[dcel.E++];
	    }

	    *em_in[0] = Edge{ein[2], &dcel.edges[2 * N - 1], em_out[0], &dcel.vertices[N], &dcel.faces[3], em_in[0] - dcel.edges};
	    *em_out[0] = Edge{&dcel.edges[N], ein[0], em_in[0], &dcel.vertices[0], &dcel.faces[2], em_out[0] - dcel.edges};
	    *em_in[1] = Edge{ein[0], &dcel.edges[N + (r + N - 1) % N], em_out[1], &dcel.vertices[N + 1], &dcel.faces[2], em_in[1] - dcel.edges};
	    *em_out[1] = Edge{&dcel.edges[N + r], ein[1], em_in[1], &dcel.vertices[r], &dcel.faces[3], em_out[1] - dcel.edges};
	    
	    ein[2]->next = dcel.edges[2 * N - 1].prev = em_in[0];
	    dcel.edges[N].next = ein[0]->prev = em_out[0];
	    ein[0]->next = dcel.edges[N + (r + N - 1) % N].prev = em_in[1];
	    dcel.edges[N + r].next = ein[1]->prev = em_out[1];
	    
        
		for(int i = 0; i < 3; ++i)
	        dcel.vertices[N + i] = Vertex{ein[i], bigTriangle[i], N + i, v_type()};
        dcel.V += 3;
        
        dcel.F = 4;
        dcel.faces[0] = Face{eout[0], 0, 0};
        dcel.faces[1] = Face{&dcel.edges[0], 1, 1};
        dcel.faces[2] = Face{ein[0], 2, 0};
        dcel.faces[3] = Face{ein[1], 3, 0};
        
        for(int i = N; i < 2 * N; ++i)
        	dcel.edges[i].adj_face = &dcel.faces[i < N + r ? 2 : 3];
        	
/*        for(int i = 0; i < dcel.E; ++i)
        	printf("%p: from %d to %d, face = %d, twin = %p\n",&dcel.edges[i], dcel.edges[i].starting_v->v_id, dcel.edges[i].next->starting_v->v_id, dcel.edges[i].adj_face->f_id, dcel.edges[i].twin);*/

        for(int i = 1; i <= 3; ++i)
            triangulate(dcel, &dcel.faces[i], i == 1);
        for(int i = 0; i < dcel.V; ++i)
        	dcel.vertices[i].v_id = i;

        for (int i = 0; i < dcel.F; ++i){
            const auto &tt = dcel.vertices_of_face(&dcel.faces[i]);
            assert((int)tt.size() == 3);
/*            for(auto p: tt)
                printf("%d ", p->v_id);
            if(dcel.faces[i].inner)
            	puts("inner");
            else
            	puts("");*/
        }

        SearchStructure ss;
        dcel.kirkpatrick_build(0, ss);

        scanf("%d", &K);
        int x, y;
        for (int i = 0; i < K; ++i) {
            scanf("%d%d", &x, &y);
            point p(x, y);
            p.rt(sc.first, sc.second);
            puts(ss.inside(p) ? "INSIDE" : "OUTSIDE");
        }
    }

    return 0;
}