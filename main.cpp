#include <iostream>
#include "read_file.h"

#include "DCEL.h"

/*
В первой строке дано целое число T — количество тестов.
Далее идут T тестов. Тесты разделены переводом строки.
N
Далее N точек — вершины многоугольника.
K
Далее K точек — запросы.
Все координаты — целые числа по модулю не превосходящие 10^9
*/

const long long INF = 1e9 + 10;
const int RAND_DENOM = 100;
const Point LEFT(-2 * INF, -INF - 1);
const Point RIGHT(2 * INF, -INF);
const Point TOP(0, 3 * INF);


void make_ccw(std::vector<Point> &points) {
    std::rotate(points.begin(), std::min_element(points.begin(), points.end()), points.end());
    const Point &a = points[1], &b = points[0], &c = points.back();
    if ((c - b) % (a - b) > 0)
        reverse(points.begin() + 1, points.end());
}


void rotate(std::vector<Point> &points, ld rot_sin, ld rot_cos) {
//    ld rot_sin = 0, rot_cos = 1;
    for (auto &p: points)
        p.rotate(rot_sin, rot_cos);
}

int main() {

    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    int T, N, K;
    scanf("%d", &T);
    ld rot_sin = (rand() % (RAND_DENOM - 1) + 1) / ld(RAND_DENOM), rot_cos = sqrtl(1 - rot_sin * rot_sin);

    std::vector<Point> big_triangle = {LEFT, RIGHT, TOP};
    rotate(big_triangle, rot_sin, rot_cos);

    for (int t = 0; t < T; ++t) {
        assert(scanf("%d", &N) == 1);
        assert(N >= 3);

        std::vector<Point> poly_points = read_points(N);

        N = (int) poly_points.size();
        make_ccw(poly_points); // the leftmost bottommost point now has index 0
        int rightmost_point_index = max_element(poly_points.begin(), poly_points.end()) - poly_points.begin();
        assert(rightmost_point_index > 0);
        rotate(poly_points, rot_sin, rot_cos);

        DCEL dcel(poly_points);
        dcel.build_triangulation_structure(rightmost_point_index, big_triangle);

        Search_Structure ss;
        dcel.kirkpatrick_build(0, ss);

        scanf("%d", &K);
        for (int i = 0; i < K; ++i) {
            int x, y;
            assert(scanf("%d%d", &x, &y) == 2);
            Point p(x, y);
            p.rotate(rot_sin, rot_cos);
            puts(ss.inside(p) ? "INSIDE" : "OUTSIDE");
        }
    }

    return 0;
}