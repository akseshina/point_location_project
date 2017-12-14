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


        std::vector<Vertex *> all_vs;
        for (int i = 0; i < my_dcel.V; ++i) {
            all_vs.push_back(&(my_dcel.vertices[i]));
        }
        my_dcel.TriangulateMonotonePolygon(all_vs);

        /*Vertex * v1 = &(my_dcel.vertices[0]);
        my_dcel.new_triangle(v1);

        Vertex * v2 = &(my_dcel.vertices[4]);
        my_dcel.new_triangle(v2);*/

        int a;

    }





    return 0;
}