#pragma once

#include <cstdio>
#include "DCEL.h"

/*
В первой строке дано целое число 𝑇 — количество тестов.
Далее идут 𝑇 тестов. Тесты разделены переводом строки.
𝑁
Далее 𝑁 точек — вершины многоугольника.
𝐾
Далее 𝐾 точек — запросы.
Все координаты — целые числа по модулю не превосходящие 10^9
*/

std::vector<point> read_points(int n) {
    std::vector<point> res;
    int x, y, cur_n;
    for (int i = 0; i < n; ++i) {
        assert(scanf("%d%d", &x, &y) == 2);
        point new_p(x, y);
        cur_n = res.size();
        // if 3 vertices are on same line
        if (cur_n >= 2 && eq(rot_matrix(res[cur_n - 2], res[cur_n - 1], new_p), 0))
            res.pop_back();
        res.push_back(point(x, y));
    }
    cur_n = res.size();
    if (eq(rot_matrix(res[cur_n - 2], res[cur_n - 1], res[0]), 0))
        res.pop_back();
    return res;
}

