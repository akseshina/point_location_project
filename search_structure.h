#pragma once

#include "geom_primitives.h"

#include <unordered_map>

struct VectorHash {
    size_t operator()(const std::vector<int> &v) const {
        size_t s = 0;
        for (int x: v)
            s = s * 239017 + (x + 1);
        return s;
    }
};

struct Search_Structure {
    std::vector<int> edges[MAX_F];
    std::vector<std::pair<Triangle, bool>> triangles;
    std::unordered_map<std::vector<int>, int, VectorHash> mapping;
    int root;

    bool inside(const Point &p) {
        return inside(p, root);
    }

    bool inside(const Point &p, int v) {
        if (edges[v].empty())
            return triangles[v].second;

        for (int u: edges[v])
            if (triangles[u].first.contains(p))
                return inside(p, u);
        return 0;
    }

    inline int get_id(const Triangle_Part &a) {
        const auto &nums = a.nums;
        int id;
        if (mapping.count(nums))
            id = mapping[nums];
        else {
            id = (int) triangles.size();
            mapping[nums] = id;
            triangles.push_back(std::make_pair(a.triangle, a.inner));
        }
        return id;
    }

    inline void add(const Triangle_Part &a, const Triangle_Part &b) {
        int a_id = get_id(a), b_id = get_id(b);
        edges[a_id].push_back(b_id);
    }
};