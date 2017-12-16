#include <iostream>
#include "DCEL.h"
#include "read_file.h"

/*
В первой строке дано целое число 𝑇 — количество тестов.
Далее идут 𝑇 тестов. Тесты разделены переводом строки.
𝑁
Далее 𝑁 точек — вершины многоугольника.
𝐾
Далее 𝐾 точек — запросы.
Все координаты — целые числа по модулю не превосходящие 10^9
*/


int main() {

    int T, N, K;
    //scanf("%d", &T);

    for (int t = 0; t < 1; ++t) { // change 1 to T
        scanf("%d", &N);
        std::vector<point> poly_points = read_points(N);
        DCEL my_dcel(poly_points);

        my_dcel.split_to_monotone(&my_dcel.faces[1]);
        //my_dcel.triangulate_monotone(&my_dcel.faces[1]);


        int a;

    }





    return 0;
}