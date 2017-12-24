#include <bits/stdc++.h>
#include "windows.h"

using namespace std;

struct Point {
    Point (int x0, int y0)
            : x(x0) ,y(y0) {}
    int x;
    int y;
};

struct Segment {
    Segment ( Point a0, Point b0)
            : a(a0) ,b(b0) {}
    Point a;
    Point b;
};

// if point q lies on line segment 'pr'
bool onSegment(Point p, Point q, Point r) {
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;

    return false;
}

// To find orientation of ordered triplet (p, q, r)
int orientation(Point p, Point q, Point r) {
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);

    if (val == 0) return 0;  // colinear

    return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool doIntersect(Point p1, Point q1, Point p2, Point q2) {
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1))
        return true;

    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1))
        return true;

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2))
        return true;

    // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2))
        return true;

    return false; // Doesn't fall in any of the above cases
}

bool doIntersect(Segment s1, Segment s2) {
    return doIntersect(s1.a, s1.b, s2.a, s2.b);
}

int rand(int a) {
    return rand() % a;
}
const int MAX = 30;
int main() {

    srand (GetTickCount());

    cout << 1<<endl;
    for(int zz = 0; zz < 1; ++zz){
    int xr, yr;

    bool done = false;
    while (!done) {

        std::vector<Segment> s;
        Point cur(rand(MAX), rand(MAX));

        s.push_back({{rand(MAX), rand(MAX)}, {rand(MAX), rand(MAX)}});
        xr = rand(MAX);
        yr = rand(MAX);
        if(s.front().a.y == s.front().b.y)
        	continue;
        map<int, bool> was;
        was[s.front().a.y] = 1;
        was[s.front().b.y] = 1;
        while (was[yr]) {
            xr = rand(MAX);
            yr = rand(MAX);
        }
        s.push_back({s.back().b, {xr, yr}});
        was[yr] = 1;

        while (true) {
            bool flag = false;

            bool attempts = true;
            int t = 0;
            while (!flag) {
                ++t;
                xr = rand(MAX);
                yr = rand(MAX);
                while (was[yr]) {
                    xr = rand(MAX);
                    yr = rand(MAX);
                }
                Segment sr(s.back().b, {xr, yr});
                flag = true;
                for (int i = 0; i < s.size() - 1; ++i) {
                    if (doIntersect(sr, s[i]))
                        flag = false;
                }
                if (flag) {
                    s.push_back(sr);
	                was[yr] = 1;
                    //cout << sr.b.x << ' ' << sr.b.y << endl;
                }

                if (t > 7) {
                    attempts = false;
                    break;
                }
            }

            if (!attempts)
                break;

            if (s.size() < 7)
                continue;

            flag = true;
            Segment send(s.back().b, s[0].a);
            if (send.a.y == send.b.y)
                continue;
            for (int i = 1; i < s.size() - 1; ++i)
                if (doIntersect(send, s[i]))
                    flag = false;
            if (flag) {
                s.push_back(send);
                cout << s.size() << endl;
                for (int i = 0; i < s.size(); ++i) {
                    cout << s[i].b.x << ' ' << s[i].b.y << endl;
                }
                cout << 676<<endl;
                for(int i = -5; i <= 20; ++i)
                for(int j = -5; j <= 20; ++j)
                	printf("%d %d\n", i, j);
/*                for (int i = 0; i < s.size(); ++i) {
                    cout << '(' << s[i].b.x << ',' << s[i].b.y << "),";
                }*/
                done = true;
                break;
            }

            if (s.size() > 15)
                break;

        }
    }
}


}