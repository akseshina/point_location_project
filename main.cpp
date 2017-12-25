#include <iostream>
#include <algorithm>

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


const int RAND_DENOM = 64;

void make_ccw(std::vector<Point> &points) {
    std::rotate(points.begin(), std::min_element(points.begin(), points.end()), points.end());
    const Point &a = points[1], &b = points[0], &c = points.back();
    if ((c - b) % (a - b) > 0)
        reverse(points.begin() + 1, points.end());
}

void rotate(std::vector<Point> &points, std::pair<ld, ld> rotate_vec) {
    for (auto &p: points)
        p.rotate(rotate_vec);
}

std::pair<ld, ld> gen_random_rotate(const std::vector<Point> &polygon, const std::vector<Point> &triangle){
    std::vector<Point> all;
    std::vector<ld> ys;
    std::pair<ld, ld> rotate_vec;
    bool ok = 0;

    int cnt = 0;
    while(!ok){
        assert(++cnt < 100);
        all = polygon;
        all.insert(all.end(), triangle.begin(), triangle.end());

        ld rot_sin = (rand() % (RAND_DENOM - 1) + 1) / ld(RAND_DENOM);
        rotate(all, rotate_vec = {rot_sin, sqrtl(1 - rot_sin * rot_sin)});

        ys.clear();
        for(const auto &p: all)
            ys.push_back(p.y);
        std::sort(ys.begin(), ys.end());

        ok = 1;
        for(int i = 0; i < int(ys.size()) - 1 && ok; ++i)
            if(eq(ys[i], ys[i + 1]))
                ok = 0;
    }
    return rotate_vec;
}

int main() {
    assert(freopen("input.txt", "r", stdin));
    assert(freopen("output.txt", "w", stdout));

    int T, N, K;
    scanf("%d", &T);

    for (int t = 0; t < T; ++t) {
        assert(scanf("%d", &N) == 1);
        assert(N >= 3);

        std::vector<Point> polygon = read_points(N);
        N = (int) polygon.size();
        std::vector<Point> polygon2 = polygon;
        
        ld INF = 0;
        for(const auto &p: polygon)
            INF = std::max(INF, std::max(std::abs(p.x), std::abs(p.y)));
        INF += 10;
        Point LEFT(-2 * INF, -INF - 1);
        Point RIGHT(2 * INF, -INF);
        Point TOP(0, 3 * INF);

        make_ccw(polygon); // the leftmost bottommost point now has index 0
        int rightmost_point_index = max_element(polygon.begin(), polygon.end()) - polygon.begin();
        assert(rightmost_point_index > 0);

        std::vector<Point> big_triangle{LEFT, RIGHT, TOP};
        std::pair<ld, ld> rotate_vec = gen_random_rotate(polygon, big_triangle);
        rotate(polygon, rotate_vec);
        rotate(big_triangle, rotate_vec);

        DCEL dcel(polygon);
        dcel.build_triangulation_structure(rightmost_point_index, big_triangle);

        Search_Structure ss;
        dcel.kirkpatrick_build(0, ss);

        scanf("%d", &K);
        for (int i = 0; i < K; ++i) {
            int x, y;
            assert(scanf("%d%d", &x, &y) == 2);
            Point p(x, y);
            p.rotate(rotate_vec);
            puts(ss.inside(p) ? "INSIDE" : "OUTSIDE");
        }
    }

    return 0;
}