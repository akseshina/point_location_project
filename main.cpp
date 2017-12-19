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


void make_ccw(std::vector<point> &points){
	std::rotate(points.begin(), std::min_element(points.begin(), points.end()), points.end());
	const point &a = points[1], &b = points[0], &c = points.back();
	if((c - b) % (a - b) > 0)
		reverse(points.begin() + 1, points.end());
}

const long long INF = 1e9 + 10;
const point LEFT(-2 * INF, -INF), RIGHT(2 * INF, -INF), TOP(0, 3 * INF);

std::vector<point> outer(const std::vector<point> &points){
	std::vector<point> p;
	p.push_back(LEFT);
	p.push_back(points[0]);
	for(int i = points.size() - 1; i >= 0; --i)
		p.push_back(points[i]);
	p.push_back(LEFT);
	--p.back().y;
	p.push_back(RIGHT);
	p.push_back(TOP);
	return p;
}

void triangulate(DCEL &dcel){
    dcel.split_to_monotone(&dcel.faces[1]);
	puts("split is ok"),fflush(stdout);
    int old_F = dcel.F;
    for (int i = 1; i < old_F; ++i)
        dcel.triangulate_monotone(&dcel.faces[i]);
}

std::pair<ld, ld> rt(std::vector<point> &points){
	ld s = rand() % 100 / 100.0, c = sqrtl(1 - s * s);
	for(auto &p: points)
		p.rt(s, c);
	return {s, c};
}

int main() {
	//freopen("input.txt", "r", stdin);

    int T, N, K;
    // scanf("%d", &T);

    for (int t = 0; t < 1; ++t) { // change 1 to T
        scanf("%d", &N);
        assert(N >= 3);
        std::vector<point> poly_points = read_points(N);
        auto sc = rt(poly_points);
/*		for(auto &p: poly_points)
			printf("%.3f %.3f\n", (double)p.x, (double)p.y);*/
		
        make_ccw(poly_points); // leftest bottoms point now has index 0
        DCEL inner_dcel(poly_points), outer_dcel(outer(poly_points));
        Edge *ein[3] = {&outer_dcel.edges[N + 4], &outer_dcel.edges[N + 3], &outer_dcel.edges[N + 2]};

        triangulate(inner_dcel);
        puts("test"), fflush(stdout);
		triangulate(outer_dcel);
        puts("test"), fflush(stdout);		
		
		std::unordered_map<const Vertex*, Vertex*> vertex_mapping;
		std::unordered_map<const Edge*, Edge*> edge_mapping;
		std::unordered_map<const Face*, Face*> face_mapping;

		DCEL dcel;
		dcel.V = N + 3;
		for(int i = 0; i < N; ++i){
			auto &v = dcel.vertices[i];
			const auto &vold1 = inner_dcel.vertices[i], &vold2 = outer_dcel.vertices[N - i];
			assert(vold1.coord == poly_points[i]);
			assert(vold2.coord == poly_points[i]);
			v.coord = poly_points[i];
			v.v_id = i;
			v.one_starting_e = vold1.one_starting_e;
			vertex_mapping[&vold1] = &v;
			vertex_mapping[&vold2] = &v;
		}
		for(int i = 0; i < 3; ++i){
			auto &v = dcel.vertices[i];
			const auto &vold = outer_dcel.vertices[N + 2 + i];
			v.coord = vold.coord;
			v.v_id = i;
			v.one_starting_e = vold.one_starting_e;
			vertex_mapping[&vold] = &v;
		}
		vertex_mapping[&outer_dcel.vertices[0]] = &dcel.vertices[N];
		vertex_mapping[&outer_dcel.vertices[N + 1]] = &dcel.vertices[0];
		
		for(int i = 0; i < inner_dcel.E; ++i){
			if(N <= i && i < 2 * N) // skip outer cycle of inner triangulation
				continue;
			auto &e = dcel.edges[dcel.E];
			const auto &eold = inner_dcel.edges[i];
			e.e_id = dcel.E;
			e.adj_face = eold.adj_face, e.prev = eold.prev, e.next = eold.next, e.twin = eold.twin, e.starting_v = eold.starting_v;
			edge_mapping[&eold] = &e;
			++dcel.E;
		}
		int M = N + 5;
		for(int i = 0; i < outer_dcel.E; ++i){
			if(M <= i && i < 2 * M) // skip outer cycle of outer triangulation (i. e. interior of original polygon + outer face)
				continue;
			auto &e = dcel.edges[dcel.E];
			const auto &eold = outer_dcel.edges[i];
			e.e_id = dcel.E;
			e.adj_face = eold.adj_face, e.prev = eold.prev, e.next = eold.next, e.twin = eold.twin, e.starting_v = eold.starting_v;
			edge_mapping[&eold] = &e;
			++dcel.E;
		}
		
		Edge *eout[3];
		for(int i = 0; i < 3; ++i){
			ein[i] = edge_mapping[ein[i]];
			assert(ein[i] != NULL);
			eout[i] = &dcel.edges[dcel.E];
			++dcel.E;
		}

		for(int i = 0; i < 3; ++i){
			auto &e = *eout[i];
			e.e_id = dcel.E;
			e.prev = eout[(i + 2) % 3], e.next = eout[(i + 1) % 3];
			e.adj_face = 0, e.twin = ein[i], e.starting_v = e.twin->next->starting_v;
		}
		

		{
			auto &f = dcel.faces[dcel.F];
			f.one_border_e = eout[0];
			f.f_id = dcel.F;
			f.inner = 0;
			++dcel.F;
		}

		for(int i = 1; i < inner_dcel.F; ++i){
			auto &f = dcel.faces[dcel.F];
			const auto &fold = inner_dcel.faces[i];
			f.one_border_e = fold.one_border_e;
			f.f_id = dcel.F;
			f.inner = 1;
			face_mapping[&fold] = &f;
			++dcel.F;
		}
		
		for(int i = 1; i < outer_dcel.F; ++i){
			auto &f = dcel.faces[dcel.F];
			const auto &fold = outer_dcel.faces[i];
			f.one_border_e = fold.one_border_e;
			f.f_id = dcel.F;
			f.inner = 0;
			face_mapping[&fold] = &f;
			++dcel.F;
		}
		
		for(int i = 0; i < dcel.V; ++i){
			auto &v = dcel.vertices[i];
			v.one_starting_e = edge_mapping[v.one_starting_e];
			assert(v.one_starting_e != NULL);
		}
		for(int i = 0; i < dcel.E; ++i){
			auto &e = dcel.edges[i];
			e.adj_face = face_mapping[e.adj_face];
			e.prev = edge_mapping[e.prev], e.next = edge_mapping[e.prev], e.twin = edge_mapping[e.twin];
			e.starting_v = vertex_mapping[e.starting_v];
			assert(e.starting_v != NULL);
		}
		for(int i = 0; i < dcel.F; ++i){
			auto &f = dcel.faces[i];
			f.one_border_e = edge_mapping[f.one_border_e];
			assert(f.one_border_e != NULL);
		}

		for(int i = 0; i < dcel.F; ++i)
			assert((int)dcel.vertices_of_face(&dcel.faces[i]).size() == 3);

        SearchStructure ss;
        dcel.kirkpatrick_build(0, ss);
        
        scanf("%d", &K);
        int x, y;
        for(int i = 0; i < K; ++i){
        	scanf("%d%d", &x, &y);
        	point p(x, y);
        	p.rt(sc.first, sc.second);
        	puts(ss.inside(p) ? "INSIDE" : "OUTSIDE");
        }
    }

    return 0;
}