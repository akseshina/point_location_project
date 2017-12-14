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
    int x, y;
    for (int i = 0; i < n; ++i) {
        scanf("%d %d", &x, &y);
        res.push_back(point(x, y));
    }
    return res;
}

