//
// Created by andrea on 07/12/23.
//

#ifndef TRUNK_LIBRARIES_H
#define TRUNK_LIBRARIES_H
#include <stdint.h>
#include <iostream>
#include <random>
#include <thread>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <utility>
#include "progressBar.h"
#include "networkit/graph/Graph.hpp"
#include "networkit/graph/Graph.hpp"
//#include "networkit/distance/Dijkstra.hpp"
#include <string>
#include "mytimer.h"
#include "networkit/centrality/DegreeCentrality.hpp"
#include "networkit/centrality/EstimateBetweenness.hpp"
#include "networkit/centrality/ApproxCloseness.hpp"
using vertex = int;
using dist = uint32_t;
const vertex null_vertex = round(std::numeric_limits<vertex>::max()/2);
const dist null_distance = round(std::numeric_limits<dist>::max()/2);

struct PQComparator
{
    bool operator() (const std::pair<dist, vertex>& p1, const std::pair<dist, vertex>& p2)
    {
        return p1.first > p2.first;
    }
};

struct PQKeywordsComparator
{
    bool operator() (const std::pair<dist, std::pair<vertex, bool>>& p1, const std::pair<dist, std::pair<vertex, bool>>& p2)
    {
        return p1.first > p2.first;
    }
};

struct PQKeywords
{
    vertex v;
    dist d1;
    dist d2;

    PQKeywords(vertex u, dist d, dist dd) : v(u), d1(d), d2(dd)
    {
    }

    bool operator<(const struct PQKeywords& other) const
    {
        //Your priority logic goes here
        return d1 > other.d1;
    }
};


template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

#endif //TRUNK_LIBRARIES_H
