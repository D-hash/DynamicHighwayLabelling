#ifndef HGHWAY_LABELING_H_
#define HGHWAY_LABELING_H_

#include <map>
#include "libraries.h"

class HighwayLabelling {
 public:
  HighwayLabelling(NetworKit::Graph &g, int l, int ordering_type, int changes, int dyntype);
    void ConstructUnweightedHighwayLabelling();
    void ConstructWeightedHighwayLabelling();
    void ConstructWeightedHighwayLabellingWithFlag();
    void ConstructDirWeighHL();
    //void RecursivePathFlag(NetworKit::Dijkstra* sssp, vertex source, vertex current_vertex, uint8_t* pruning_flag);
    void ConstructUnweightedNaiveLabelling();
    void ConstructWeightedNaiveLabelling();
    void ConstructDirWeighNL();
  int GetNumberOfNodes(){return V;};
  long LabellingSize();
  long DirectedLabellingSize();
  long NaiveLabellingSize();
  long DirectedNaiveLabellingSize();
  dist min(dist a, dist b);

  // Returns distance vetween vertices v and w if they are connected.
  dist QueryDistance(vertex s, vertex t);
    dist DirectedQueryDistance(vertex s, vertex t);
  dist QueryDistanceQuadratic(vertex s, vertex t);
  dist NaiveQueryDistance(vertex s, vertex t);
  dist DirectedNaiveQueryDistance(vertex s, vertex t);
  dist BFSQuery(vertex s, vertex t);
    dist DijkstraQuery(vertex s, vertex t);
    dist DirectedDijkstra(vertex s, vertex t);
  dist BoundedSearch(vertex s, vertex t);
  dist P2PBFS(vertex s, vertex landmark, dist ub, dist partial, dist* to_vertices);
    void GetBeerStores(std::vector<vertex>& beer_stores);
    void GetIncrementalBeerStores(std::vector<vertex>& beer_stores);
  void SetBeerStores(std::vector<vertex>& beer_stores);
  void StoreIndex(std::string filename, bool dynamic);

  // INCREMENTAL
  void AddLandmark(vertex r);
    void AddLandmarkDirected(vertex r);

  // DECREMENTAL
  void RemoveLandmark(vertex r);
    void RemoveLandmarkDirected(vertex r);
  std::vector<vertex> reverse_ordering;
  vertex L;
  std::set<vertex> landmarks;

private:
  vertex V;  // total number of vertices
  vertex E; // total number of edges
   // total number of landmarks
  NetworKit::Graph &graph;
    std::vector<std::vector<vertex>> landmarks_distances;
    std::vector<std::vector<vertex>> in_landmarks_distances;
    std::vector<std::vector<vertex>> out_landmarks_distances;
    std::vector<std::vector<dist>> distances;
    std::vector<std::vector<dist>> in_distances;
    std::vector<std::vector<dist>> out_distances;
  std::unordered_map<vertex, std::unordered_map<vertex, dist>> highway;
  std::vector<vertex> ordering;
  std::vector<std::vector<dist>> naive_labeling;
  std::vector<std::vector<dist>> in_naive_labeling;
  std::vector<std::vector<dist>> out_naive_labeling;
  std::vector<vertex> landmarks_incremental;
    // temp structures
    std::vector<bool> pruning_flag;
    std::vector<bool> settled;
    std::vector<dist> dij_distances;
    std::vector<dist> dij_distances_in;
    std::vector<vertex> predecessors;
};


HighwayLabelling::HighwayLabelling(NetworKit::Graph &g, int l, int ordering_type, int changes, int dyntype)
    : graph(g) {
    V = 0;
    E = 0;
    L = l;
    graph = g;

    V = this->graph.numberOfNodes();
    
    E = this->graph.numberOfEdges();
    auto *ordering_rank = new std::pair<double, vertex>[graph.numberOfNodes()];
    this->ordering.resize(graph.numberOfNodes());
    this->reverse_ordering.resize(graph.numberOfNodes());

    this->graph.parallelForNodes([&](vertex i) {
        assert(graph.hasNode(i));
        this->ordering[i] = null_vertex;
        this->reverse_ordering[i] = null_vertex;
        ordering_rank[i] = {0, i};
    });


    double centr_time = 0.0;

    if (ordering_type == 0) {
        INFO("BY DEGREE");
        NetworKit::DegreeCentrality *rank = new NetworKit::DegreeCentrality(graph);
        rank->run();
        this->graph.forNodes([&](vertex i) {
            assert(graph.hasNode(i));
            ordering_rank[i] = std::make_pair(rank->score(i), i);
        });
        delete rank;
    }

    if (ordering_type == 1) {
        INFO("RANDOM ORDERING");
        std::vector<vertex> t(graph.numberOfNodes());
        for (vertex i = 0; i < graph.numberOfNodes(); i++) t[i] = i;
        std::shuffle(t.begin(), t.end(), std::mt19937());
        for (vertex i = 0; i < graph.numberOfNodes(); i++) ordering_rank[i] = std::make_pair(t[i], i);
    }

    if (ordering_type == 2) {
        INFO("BY APX BETW");
        mytimer local_constr_timer;
        double max_time = 30.0;
        double cumulative_time = 0.0;
        double fract = 0.33;
        double n_samples = round(std::pow((double) graph.numberOfNodes(), fract));

        while (cumulative_time < max_time && n_samples < (double) graph.numberOfNodes()) {
            local_constr_timer.restart();

            std::cout << "fract: " << fract << " " << n_samples << " SAMPLES\n";
            NetworKit::EstimateBetweenness *rank = new NetworKit::EstimateBetweenness(graph, n_samples, false, true);

            rank->run();

            this->graph.forNodes([&](vertex i) {
                assert(graph.hasNode(i));
                assert(i<graph.numberOfNodes());
                ordering_rank[i] = std::make_pair(rank->score(i), i);
            });
            delete rank;
            cumulative_time += local_constr_timer.elapsed();
            n_samples *= 2;
        }
    }
    if (ordering_type == 3) {
        srandom(2024);
        INFO("BY DISTANCE-" + std::to_string((vertex) log2(V)) +" BOUNDED DOMINATING SET");
        std::set<vertex> to_cover;
        for (vertex v: graph.nodeRange()) to_cover.insert(v);
        vertex landmark_counter = 0;
        while (!to_cover.empty() && landmark_counter < L) {
            std::queue<std::pair<vertex, dist> > que;
            auto it = to_cover.begin();
            std::advance(it, random() % to_cover.size());
            vertex source = *it;
            landmark_counter++;
            que.push(std::make_pair(source, 0));
            while (!que.empty()) {
                auto p = que.front();
                que.pop();
                if (p.second > ((vertex) log2(V))) break;
                ordering_rank[p.first] = std::make_pair(to_cover.size(), p.first);
                to_cover.erase(p.first);
                for (vertex v: graph.neighborRange(p.first)) {
                    if (to_cover.find(v) != to_cover.end())
                        que.push(std::make_pair(v, p.second + 1));
                }
            }
            //ordering_rank[source] = std::make_pair(graph.numberOfNodes()+1+landmark_counter, source);
        }
    }

    if (ordering_type == 4) {
        INFO("BY APPROX CLOSENESS");
        NetworKit::ApproxCloseness *rank = new NetworKit::ApproxCloseness(graph, 1000);
        rank->run();
        this->graph.forNodes([&](vertex i) {
            assert(graph.hasNode(i));
            ordering_rank[i] = std::make_pair(rank->score(i), i);
        });
        delete rank;
    }

    std::sort(ordering_rank, ordering_rank + graph.numberOfNodes(),
              [](const std::pair<double, vertex> &a, const std::pair<double, vertex> &b) {
                  if (a.first == b.first)
                      return a.second > b.second;
                  else {
                      return a.first > b.first;
                  }
              });


    for (size_t count = 0; count < graph.numberOfNodes(); count++) {
        this->reverse_ordering[count] = ordering_rank[count].second;
        this->ordering[ordering_rank[count].second] = count;
    }
    landmarks.clear();
    if(dyntype == 1) {
        for (vertex i = 0; i < L - changes ; i++) {
            landmarks.insert(reverse_ordering[i]);
        }
        for(vertex i = L - changes; i < L; i++){
            landmarks_incremental.push_back(reverse_ordering[i]);
        }
        L = L - changes;
    }
    else if(dyntype == 2) {
        for (vertex i = 0; i < L ; i++) {
            landmarks.insert(reverse_ordering[i]);
        }
    }
    delete[] ordering_rank;

    pruning_flag.resize(V);
    dij_distances.resize(V, null_distance);
    dij_distances_in.resize(V, null_distance);
    predecessors.resize(V, null_vertex);
    settled.resize(V, false);
}

void HighwayLabelling::GetBeerStores(std::vector<vertex> &beer_stores) {
    for(auto v: landmarks){
        beer_stores.push_back(v);
    }
}

void HighwayLabelling::GetIncrementalBeerStores(std::vector<vertex> &beer_stores) {
    for(auto v: landmarks_incremental){
        beer_stores.push_back(v);
    }
}

void HighwayLabelling::SetBeerStores(std::vector<vertex> &beer_stores) {
    landmarks.clear();
    L = beer_stores.size();
    for(auto v: beer_stores){
        landmarks.insert(v);
    }
}

long HighwayLabelling::LabellingSize() {
  long size = 0;
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < distances[i].size(); j++) {
      if(distances[i][j] != null_distance)
        size++;
    }
  }

  for (const vertex & v: landmarks) {
    for (const vertex & w: landmarks) {
      if(highway[v][w] != null_distance)
       size++;
    }
  }

  return size;
}

long HighwayLabelling::DirectedLabellingSize() {
  long size = 0;
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < in_distances[i].size(); j++) {
      if(in_distances[i][j] != null_distance)
        size++;
    }
    for (int j = 0; j < out_distances[i].size(); j++) {
      if(out_distances[i][j] != null_distance)
        size++;
    }
  }

  for (const vertex & v: landmarks) {
    for (const vertex & w: landmarks) {
      if(highway[v][w] != null_distance)
       size++;
    }
  }

  return size;
}

long HighwayLabelling::NaiveLabellingSize() {
    long size = 0;
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < L; j++) {
            if(naive_labeling[i][j] != null_distance)
                size++;
        }
    }

    return size;
}

long HighwayLabelling::DirectedNaiveLabellingSize() {
    long size = 0;
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < L; j++) {
            if(in_naive_labeling[i][j] != null_distance)
                size++;
            if(out_naive_labeling[i][j] != null_distance)
                size++;
        }
    }

    return size;
}

void HighwayLabelling::ConstructUnweightedHighwayLabelling() {}
//void HighwayLabelling::ConstructUnweightedHighwayLabelling() {
//  // Initialization
//    distances.resize(V);
//  for(int i = 0; i < V; i++) {
//      distances[i].resize(L);
//    for(int j = 0; j < L; j++)
//      distances[i][j] = null_distance;
//  }
//
//  highway.resize(L);
//  for(int i = 0; i < L; i++)
//    highway[i].resize(L);
//
//  // Start computing Highway Labelling (HL)
//  ProgressStream hl_bar(L);
//
//  hl_bar.label() << "Unweighted highway labeling construction";
//  for (int i = 0; i < L; i++) {
//    dist *P = new dist[V];
//    for(int j = 0; j < V; j++)
//      P[j] = null_distance;
//
//    std::queue<int> que[2];
//
//    que[0].push(reverse_ordering[i]); que[0].push(-1);
//    distances[reverse_ordering[i]][i] = 0; P[reverse_ordering[i]] = 0; int use = 0;
//    while (!que[0].empty()) {
//      int u = que[use].front();
//      que[use].pop();
//
//      if(u == -1) {
//        use = 1 - use;
//        que[use].push(-1);
//        continue;
//      }
//
//      for (vertex w : graph.neighborRange(u)) {
//        if (P[w] == null_distance) {
//          P[w] = P[u] + 1;
//          if(use == 1 || ordering[w] < L)
//            que[1].push(w);
//          else {
//            que[0].push(w);
//            distances[w][i] = P[w];
//          }
//        }
//      }
//    }
//
//    for(int j = 0; j < L; j++) {
//      if(P[reverse_ordering[j]] != null_distance) {
//        highway[i][j] = P[reverse_ordering[j]];
//        highway[j][i] = P[reverse_ordering[j]];
//      }
//    }
//
//    delete [] P;
//    ++hl_bar;
//  }
//}



void HighwayLabelling::ConstructWeightedHighwayLabellingWithFlag() {
    // Initialization
    distances.resize(V);
    landmarks_distances.resize(V);
    for(int i = 0; i < V; i++) {
        distances[i].clear();
        landmarks_distances[i].clear();
    }
    highway.clear();
    for(const auto & l1: landmarks) {
        highway[l1] = std::unordered_map<vertex, dist> ();
    }

    // Start computing Highway Labelling (HL)
    ProgressStream hl_bar(L);

    hl_bar.label() << "Weighted highway labeling construction";
    std::vector<vertex> reached_vertices;
    for(const vertex & b : landmarks) {

        highway[b][b] = 0;
        dij_distances[b] = 0;
        reached_vertices.push_back(b);
        landmarks_distances[b].push_back(b);
        distances[b].push_back(0);
        std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
                PQComparator> pq;
        pq.push(std::make_pair(0,b));
        while(!pq.empty()){
            vertex v = pq.top().second;
            pq.pop();
            if(landmarks.find(v) != landmarks.end() && v != b){
                pruning_flag[v] = true;
                highway[b][v] = dij_distances[v];
                highway[v][b] = dij_distances[v];
            }
            if(!pruning_flag[v] && v != b){
                std::vector<vertex>::iterator insertion_index;
                for(size_t i = 0; i < landmarks_distances[v].size(); i++) {
                    if (landmarks_distances[v][i] == b) {
                        goto check;
                    }
                }
                insertion_index = std::upper_bound(landmarks_distances[v].begin(),
                                                        landmarks_distances[v].end(), b);
                distances[v].insert(
                        distances[v].begin() + (insertion_index - landmarks_distances[v].begin()),
                        dij_distances[v]);
                landmarks_distances[v].insert(insertion_index, b);
                check:{};
            }
            for(auto w: graph.neighborRange(v)){
                if(dij_distances[w] > dij_distances[v] + graph.weight(v,w)){
                    dij_distances[w] = dij_distances[v] + (dist)graph.weight(v,w);
                    predecessors[w] = v;
                    pruning_flag[w] = pruning_flag[v];
                    reached_vertices.push_back(w);
                    pq.push(std::make_pair(dij_distances[w], w));
                }

            }
        }
        for(const auto & v: reached_vertices){
            dij_distances[v] = null_distance;
            predecessors[v] = null_vertex;
            pruning_flag[v] = false;
        }
        reached_vertices.clear();
        ++hl_bar;
    }

}

void HighwayLabelling::ConstructDirWeighHL() {
    // Initialization
    in_distances.resize(V);
    out_distances.resize(V);
    in_landmarks_distances.resize(V);
    out_landmarks_distances.resize(V);
    for(int i = 0; i < V; i++) {
        in_distances[i].clear();
        out_distances[i].clear();
        in_landmarks_distances[i].clear();
        out_landmarks_distances[i].clear();
    }
    highway.clear();
    for(const auto & l1: landmarks) {
        highway[l1] = std::unordered_map<vertex, dist> ();
        for(const auto & l2: landmarks) {
            highway[l1][l2] = null_distance;
        }
    }

    // Start computing Highway Labelling (HL)
    ProgressStream hl_bar(L);

    hl_bar.label() << "Directed weighted highway labeling construction";
    std::vector<vertex> reached_vertices;
    // FORWARD
    for(const vertex & b : landmarks) {
        vertex lndm_to_reach = landmarks.size()-1;
        highway[b][b] = 0;
        dij_distances[b] = 0;
        reached_vertices.push_back(b);
        out_landmarks_distances[b].push_back(b);
        out_distances[b].push_back(0);
        std::priority_queue<std::pair<dist, vertex>, std::vector<std::pair<dist, vertex>>,
                PQComparator> pq;
        pq.push(std::make_pair(0, b));
        while (!pq.empty()) {
            vertex v = pq.top().second;
            pq.pop();
            if(settled[v]) continue;
            if(lndm_to_reach == 0 && DirectedQueryDistance(b, v) <= dij_distances[v]) continue;
            if (landmarks.find(v) != landmarks.end() && v != b) {
                pruning_flag[v] = true;
                highway[b][v] = dij_distances[v];
                lndm_to_reach--;
                settled[v] = true;
            }
            if (!pruning_flag[v] && v != b && !settled[v]) {
                settled[v] = true;
                out_distances[v].push_back(dij_distances[v]);
                out_landmarks_distances[v].push_back(b);
            }
            for (auto w: graph.neighborRange(v)) {
                if (dij_distances[w] > dij_distances[v] + graph.weight(v, w)) {
                    dij_distances[w] = dij_distances[v] + static_cast<dist>(graph.weight(v, w));
                    predecessors[w] = v;
                    pruning_flag[w] = pruning_flag[v];
                    reached_vertices.push_back(w);
                    pq.push(std::make_pair(dij_distances[w], w));
                }

            }
        }
        for (const auto &v: reached_vertices) {
            dij_distances[v] = null_distance;
            predecessors[v] = null_vertex;
            pruning_flag[v] = false;
            settled[v] = false;
        }
        reached_vertices.clear();
        lndm_to_reach = landmarks.size() - 1;
        // REVERSE
        dij_distances[b] = 0;
        reached_vertices.push_back(b);
        in_landmarks_distances[b].push_back(b);
        in_distances[b].push_back(0);
        while(!pq.empty()) pq.pop();
        pq.push(std::make_pair(0, b));
        while (!pq.empty()) {
            vertex v = pq.top().second;
            pq.pop();
            if(settled[v]) continue;
            if(lndm_to_reach == 0 && DirectedQueryDistance(v, b) <= dij_distances[v]) continue;
            if (landmarks.find(v) != landmarks.end() && v != b) {
                pruning_flag[v] = true;
                highway[v][b] = dij_distances[v];
                lndm_to_reach--;
                settled[v] = true;
            }
            if (!pruning_flag[v] && v != b && !settled[v]) {
                in_distances[v].push_back(dij_distances[v]);
                settled[v] = true;
                in_landmarks_distances[v].push_back(b);
            }
            for (auto w: graph.inNeighborRange(v)) {
                if (dij_distances[w] > dij_distances[v] + graph.weight(w,v)) {
                    dij_distances[w] = dij_distances[v] + (dist) graph.weight(w,v);
                    predecessors[w] = v;
                    pruning_flag[w] = pruning_flag[v];
                    reached_vertices.push_back(w);
                    pq.push(std::make_pair(dij_distances[w], w));
                }

            }
        }
        for (const auto &v: reached_vertices) {
            dij_distances[v] = null_distance;
            predecessors[v] = null_vertex;
            pruning_flag[v] = false;
            settled[v] = false;
        }
        reached_vertices.clear();
        ++hl_bar;
    }
}

void HighwayLabelling::ConstructUnweightedNaiveLabelling() {
    // Initialization
    naive_labeling.resize(V);
    for(vertex i = 0; i < V; i++) {
        naive_labeling[i].resize(L);
        for(vertex j = 0; j < L; j++)
            naive_labeling[i][j] = null_distance;
    }


    // Start computing Naive Labelling
    ProgressStream na_bar(L);

    na_bar.label() << "Unweighted naive labeling construction";
    for(vertex s = 0; s < L; s++){
        dist *P = new dist[V];
        for(vertex j = 0; j < V; j++)
            P[j] = null_distance;
        std::queue<vertex> que;
        que.push(reverse_ordering[s]);
        naive_labeling[reverse_ordering[s]][s] = 0; P[reverse_ordering[s]] = 0;

        while(!que.empty()){
            vertex v = que.front();
            que.pop();
            naive_labeling[v][s] = P[v];

            for(vertex w: graph.neighborRange(v)){
                if(P[w]==null_distance){
                    P[w] = P[v] + 1;
                    que.push(w);
                }
            }
        }
    ++na_bar;
    }

}

void HighwayLabelling::ConstructWeightedNaiveLabelling() {
    // Initialization
    naive_labeling.resize(V);
    for(vertex i = 0; i < V; i++) {
        naive_labeling[i].resize(L);
        for(vertex j = 0; j < L; j++)
            naive_labeling[i][j] = null_distance;
    }
    ProgressStream na_bar(L);

    na_bar.label() << "Weighted naive labeling construction";
    for(vertex b = 0; b < L; b++) {
        bool* pruning_flag = new bool[V];

        naive_labeling[reverse_ordering[b]][b] = 0;
        std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
                PQComparator> pq;
        pq.push(std::make_pair(0,reverse_ordering[b]));
        while(!pq.empty()){
            vertex v = pq.top().second;
            pq.pop();

            for(auto w: graph.neighborRange(v)){
                if(naive_labeling[w][b] > naive_labeling[v][b] + graph.weight(v,w)){
                    naive_labeling[w][b] = naive_labeling[v][b] + (dist)graph.weight(v,w);

                    pruning_flag[w] = pruning_flag[v];
                    pq.push(std::make_pair(naive_labeling[w][b], w));
                }

            }
        }
        ++na_bar;
        delete [] pruning_flag;
    }
    
}

void HighwayLabelling::ConstructDirWeighNL() {
    // Initialization
    in_naive_labeling.resize(V);
    out_naive_labeling.resize(V);
    for(vertex i = 0; i < V; i++) {
        in_naive_labeling[i].resize(L);
        out_naive_labeling[i].resize(L);
        for(vertex j = 0; j < L; j++){
            in_naive_labeling[i][j] = null_distance;
            out_naive_labeling[i][j] = null_distance;
            }
    }
    ProgressStream na_bar(L);

    na_bar.label() << "Directed weighted naive labeling construction";
    for(vertex b = 0; b < L; b++) {
        bool* pruning_flag = new bool[V];
        std::vector<vertex> reached_vertices;
        out_naive_labeling[reverse_ordering[b]][b] = 0;
        std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
                PQComparator> pq;
        pq.push(std::make_pair(0,reverse_ordering[b]));
        while(!pq.empty()){
            vertex v = pq.top().second;
            pq.pop();
            reached_vertices.push_back(v);
            for(auto w: graph.neighborRange(v)){
                if(out_naive_labeling[w][b] > out_naive_labeling[v][b] + graph.weight(v,w)){
                    out_naive_labeling[w][b] = out_naive_labeling[v][b] + (dist)graph.weight(v,w);

                    pruning_flag[w] = pruning_flag[v];
                    pq.push(std::make_pair(out_naive_labeling[w][b], w));
                }

            }
        }

        for(const auto & v: reached_vertices){
            pruning_flag[v] = false;
        }
        reached_vertices.clear();
        in_naive_labeling[reverse_ordering[b]][b] = 0;
        while(!pq.empty()) pq.pop();
        pq.push(std::make_pair(0,reverse_ordering[b]));
        while(!pq.empty()){
            vertex v = pq.top().second;
            pq.pop();

            for(auto w: graph.inNeighborRange(v)){
                if(in_naive_labeling[w][b] > in_naive_labeling[v][b] + graph.weight(w,v)){
                    in_naive_labeling[w][b] = in_naive_labeling[v][b] + (dist)graph.weight(w,v);

                    pruning_flag[w] = pruning_flag[v];
                    pq.push(std::make_pair(in_naive_labeling[w][b], w));
                }

            }
        }

        ++na_bar;
        delete [] pruning_flag;
    }
    
}

dist HighwayLabelling::min(dist a, dist b) {
  return (a < b) ? a : b;
}

dist HighwayLabelling::QueryDistance(vertex s, vertex t) {

  dist m = null_distance, uni1[landmarks_distances[s].size()] = {0}, uni2[landmarks_distances[t].size()] = {0}; vertex i = 0, j = 0, i1 = 0, j1 = 0;
  while (i < landmarks_distances[s].size() && j < landmarks_distances[t].size()) {
    if (landmarks_distances[s][i] < landmarks_distances[t][j]) {
      uni1[i1] = i; i++; i1++;
    } else if (landmarks_distances[t][j] < landmarks_distances[s][i]) {
      uni2[j1] = j; j++; j1++;
    } else {
      m = min(m, distances[s][i] + distances[t][j]);
      i++; j++;
    }
  }

  while (i < landmarks_distances[s].size()) {
    uni1[i1] = i; i++; i1++;
  }

  while (j < landmarks_distances[t].size()) {
    uni2[j1] = j; j++; j1++;
  }

  i = 0;
  while (i < i1) { j = 0;
    while (j < j1) {
      m = min(m, distances[s][uni1[i]] + highway[landmarks_distances[s][uni1[i]]][landmarks_distances[t][uni2[j]]] + distances[t][uni2[j]]);
      j++;
    }
    i++;
  }

  return m;
}

dist HighwayLabelling::DirectedQueryDistance(vertex s, vertex t) {

    dist m = null_distance, uni1[in_landmarks_distances[s].size()] = {0}, uni2[out_landmarks_distances[t].size()] = {0}; vertex i = 0, j = 0, i1 = 0, j1 = 0;
    while (i < in_landmarks_distances[s].size() && j < out_landmarks_distances[t].size()) {
        if (in_landmarks_distances[s][i] < out_landmarks_distances[t][j]) {
            uni1[i1] = i; i++; i1++;
        } else if (out_landmarks_distances[t][j] < in_landmarks_distances[s][i]) {
            uni2[j1] = j; j++; j1++;
        } else {
            m = min(m, in_distances[s][i] + out_distances[t][j]);
            i++; j++;
        }
    }

    while (i < in_landmarks_distances[s].size()) {
        uni1[i1] = i; i++; i1++;
    }

    while (j < out_landmarks_distances[t].size()) {
        uni2[j1] = j; j++; j1++;
    }

    i = 0;
    // while (i < i1) { j = 0;
    //     while (j < j1) {
    //         m = min(m, in_distances[s][uni1[i]] + highway[in_landmarks_distances[s][uni1[i]]][out_landmarks_distances[t][uni2[j]]] + out_distances[t][uni2[j]]);
    //         j++;
    //     }
    //     i++;
    // }
    while (i < i1) {
        j = 0;
        for(j = 0; j < out_landmarks_distances[t].size(); j++) {
            m = min(m, in_distances[s][uni1[i]] +
                highway[in_landmarks_distances[s][uni1[i]]][out_landmarks_distances[t][j]] +
                out_distances[t][j]);
        }
        i++;
    }

    j = 0;
    while (j < j1) {
        i = 0;
        for(i = 0; i < in_landmarks_distances[s].size(); i++) {
            m = min(m, in_distances[s][i] +
                highway[in_landmarks_distances[s][i]][out_landmarks_distances[t][uni2[j]]] +
                out_distances[t][uni2[j]]);
        }
        j++;
    }

    return m;
}

dist HighwayLabelling::QueryDistanceQuadratic(vertex s, vertex t) {
    dist m = null_distance;
    vertex i, j;
    for(i = 0; i < L; i++) {
        for (j = 0; j < L; j++)
            m = min(m, distances[s][i] + highway[landmarks_distances[s][i]][landmarks_distances[t][j]] + distances[t][j]);
    }

    return m;
}

dist HighwayLabelling::NaiveQueryDistance(vertex s, vertex t) {
    dist m = null_distance;
    for(vertex i = 0; i < L; i++){
            m = min(m, naive_labeling[s][i]+naive_labeling[t][i]);
    }
    return m;
}

dist HighwayLabelling::DirectedNaiveQueryDistance(vertex s, vertex t) {
    dist m = null_distance;
    for(vertex i = 0; i < L; i++){
            m = min(m, in_naive_labeling[s][i]+out_naive_labeling[t][i]);
    }
    return m;
}

dist HighwayLabelling::BFSQuery(vertex s, vertex t) {
    // Initialization
    dist *s_to_vertices = new dist[V];
    dist *t_to_vertices = new dist[V];
    for(vertex j = 0; j < V; j++) {
        s_to_vertices[j] = null_distance;
        t_to_vertices[j] = null_distance;
    }
    std::queue<vertex> que;
    que.push(s);
    s_to_vertices[s] = 0;

    while(!que.empty()){
        vertex v = que.front();
        que.pop();

        for(vertex w: graph.neighborRange(v)){
            if(s_to_vertices[w]==null_distance){
                s_to_vertices[w] = s_to_vertices[v] + 1;
                que.push(w);
            }
        }
    }

    que.push(t);
    t_to_vertices[t] = 0;

    while(!que.empty()){
        vertex v = que.front();
        que.pop();

        for(vertex w: graph.neighborRange(v)){
            if(t_to_vertices[w]==null_distance){
                t_to_vertices[w] = t_to_vertices[v] + 1;
                que.push(w);
            }
        }
    }

    dist m = null_distance;
    for(vertex l = 0; l < L; l++){
        m = min(s_to_vertices[reverse_ordering[l]] + t_to_vertices[reverse_ordering[l]], m);
    }

    delete [] s_to_vertices;
    delete [] t_to_vertices;
    return m;
}

dist HighwayLabelling::DijkstraQuery(vertex s, vertex t) {
    dist* s_distances = new dist[V];
    for(vertex i = 0; i < V; i++){
        s_distances[i] = null_distance;
    }
    s_distances[s] = 0;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pq;
    pq.push(std::make_pair(0,s));
    while(!pq.empty()){
        vertex v = pq.top().second;
        pq.pop();
        for(auto w: graph.neighborRange(v)){
            if(s_distances[w] > s_distances[v] + graph.weight(v,w)){
                s_distances[w] = s_distances[v] + (dist)graph.weight(v,w);

                pq.push(std::make_pair(s_distances[w], w));
            }

        }
    }

    dist* t_distances = new dist[V];
    for(vertex i = 0; i < V; i++){
        t_distances[i] = null_distance;
    }
    t_distances[t] = 0;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pqt;
    pqt.push(std::make_pair(0,t));
    while(!pqt.empty()){
        vertex v = pqt.top().second;
        pqt.pop();
        for(auto w: graph.neighborRange(v)){
            if(t_distances[w] > t_distances[v] + graph.weight(v,w)){
                t_distances[w] = t_distances[v] + (dist)graph.weight(v,w);

                pqt.push(std::make_pair(t_distances[w], w));
            }

        }
    }

//    auto s_to_vertices = NetworKit::Dijkstra(graph, s, false, false);
//    s_to_vertices.run();
//    auto t_to_vertices = NetworKit::Dijkstra(graph, t, false, false);
//    t_to_vertices.run();
//    dist m = null_distance;
//    for(vertex l = 0; l < L; l++){
//        m = min(s_to_vertices.distance(reverse_ordering[l]) + t_to_vertices.distance(reverse_ordering[l]), m);
//    }
    dist m = null_distance;
    for(const vertex & l: landmarks){
        m = min(m, s_distances[l] + t_distances[l]);
    }
    delete [] s_distances;
    delete [] t_distances;
    return m;
}

dist HighwayLabelling::DirectedDijkstra(vertex s, vertex t) {
    dist* s_distances = new dist[V];
    for(vertex i = 0; i < V; i++){
        s_distances[i] = null_distance;
    }
    s_distances[s] = 0;
    std::set<vertex> rlndm;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pq;
    pq.push(std::make_pair(0,s));
    while(!pq.empty()){
        vertex v = pq.top().second;
        pq.pop();
        if(landmarks.find(v) != landmarks.end()) {
            rlndm.insert(v);
        }
        if(rlndm.size() == landmarks.size()) break;
        for(auto w: graph.neighborRange(v)){
            if(s_distances[w] > s_distances[v] + graph.weight(v,w)){
                s_distances[w] = s_distances[v] + (dist)graph.weight(v,w);

                pq.push(std::make_pair(s_distances[w], w));
            }

        }
    }

    dist* t_distances = new dist[V];
    for(vertex i = 0; i < V; i++){
        t_distances[i] = null_distance;
    }
    rlndm.clear();
    t_distances[t] = 0;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pqt;
    pqt.push(std::make_pair(0,t));
    while(!pqt.empty()){
        vertex v = pqt.top().second;
        pqt.pop();
        if(landmarks.find(v) != landmarks.end()) {
            rlndm.insert(v);
        }
        if(rlndm.size() == landmarks.size()) break;
        for(auto w: graph.inNeighborRange(v)){
            if(t_distances[w] > t_distances[v] + graph.weight(w,v)){
                t_distances[w] = t_distances[v] + (dist)graph.weight(w,v);

                pqt.push(std::make_pair(t_distances[w], w));
            }

        }
    }

    dist m = null_distance;
    vertex v = null_vertex;
    dist from_s = null_distance;
    dist to_t = null_distance;
    for(const vertex & l: landmarks){
        if(m > s_distances[l] + t_distances[l]) {
            v = l;
            m = s_distances[l] + t_distances[l];
            from_s = s_distances[l];
            to_t = t_distances[l];
        }
        //m = min(m, s_distances[l] + t_distances[l]);
    }
    delete [] s_distances;
    delete [] t_distances;
    //std::cout << "hub " << v << " from s " << from_s << " to " << to_t << "\n";
    return m;
}

dist HighwayLabelling::P2PBFS(vertex s, vertex landmark, dist ub, dist partial, dist* to_vertices){
    std::vector<vertex> visited;
    std::queue<vertex> que;
    que.push(s);
    to_vertices[s] = 0;
    visited.push_back(s);

    while(!que.empty()){
        vertex v = que.front();
        que.pop();
        if(to_vertices[v] >= ub || to_vertices[v] + partial >= ub){
            for(auto u: visited)
                to_vertices[u] = null_distance;
            return null_distance;
        }
        if(v == landmark){
            dist m = to_vertices[v];
            for(auto u: visited)
                to_vertices[u] = null_distance;
            return m + partial;
        }
        for(vertex w: graph.neighborRange(v)){
            if(to_vertices[w]==null_distance){
                to_vertices[w] = to_vertices[v] + 1;
                visited.push_back(w);

                que.push(w);
            }
        }
    }

    for(auto u: visited)
        to_vertices[u] = null_distance;
    return null_distance;
}

void HighwayLabelling::AddLandmark(vertex r) {
    // compute Highway distances from previous landmarks to r
    // add distances in Highway between r and landmarks directly connected to r (i.e. in L(r))
    landmarks.insert(r);
    highway[r] = std::unordered_map<vertex, dist> ();
    for(size_t i = 0; i < landmarks_distances[r].size(); i++){
        highway[r][landmarks_distances[r][i]] = distances[r][i];
        highway[landmarks_distances[r][i]][r] = distances[r][i];
    }
    highway[r][r] = 0;
    // fill distances in Highway between r and landmarks not covering r (i.e. not in L(r))
    for(const vertex & l: landmarks){
        if(highway[r].find(l) == highway[r].end()){
            highway[r][l] = null_distance;
            for(size_t j = 0; j < landmarks_distances[r].size(); j++){
                highway[r][l] = std::min(highway[r][l],highway[r][landmarks_distances[r][j]] + highway[landmarks_distances[r][j]][l]);
            }
            highway[l][r] = highway[r][l];
        }
    }
    landmarks_distances[r].clear();
    landmarks_distances[r].push_back(r);
    distances[r].clear();
    distances[r].push_back(0);
    // search from r to fill L a la Akiba (prune when reaching landmarks different from r, compute query to vertices not in L)
    std::vector<vertex> reached_vertices;
    dij_distances[r] = 0;
    reached_vertices.push_back(r);
    std::vector<vertex> reached_landmarks;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pq;
    for(auto w: graph.neighborRange(r)){
            dij_distances[w] = dij_distances[r] + (dist)graph.weight(r,w);
            pq.push(std::make_pair(dij_distances[w], w));
            reached_vertices.push_back(w);
        }
    while(!pq.empty()){
        vertex v = pq.top().second;
        pq.pop();
        if(pruning_flag[v]) continue;
        if(landmarks.find(v) != landmarks.end() && v != r){
            pruning_flag[v] = true;
            reached_landmarks.push_back(v);
            continue;
        }
        dist query_dist = QueryDistance(r,v);
        if(query_dist <= dij_distances[v]){
            pruning_flag[v] = true;
            continue;
        }
        auto insertion_index = std::upper_bound(landmarks_distances[v].begin(), landmarks_distances[v].end(), r);
        distances[v].insert(distances[v].begin() + (insertion_index - landmarks_distances[v].begin()), dij_distances[v]);
        landmarks_distances[v].insert(insertion_index, r);
        for(auto w: graph.neighborRange(v)){
            if(dij_distances[w] > dij_distances[v] + graph.weight(v,w)){
                dij_distances[w] = dij_distances[v] + (dist)graph.weight(v,w);
                predecessors[w] = v;
                pruning_flag[w] = pruning_flag[v];
                pq.push(std::make_pair(dij_distances[w], w));
                reached_vertices.push_back(w);
            }

        }
    }

    for(const auto & lndm: reached_landmarks){
        for(const auto & v: reached_vertices){
            for(size_t i = 0; i < landmarks_distances[v].size(); i++){
                if(landmarks_distances[v][i] < lndm) continue;
                if(landmarks_distances[v][i] > lndm) break;
                if(distances[v][i] >= dij_distances[v] + highway[r][lndm]){
                    landmarks_distances[v].erase(landmarks_distances[v].begin()+i);
                    distances[v].erase(distances[v].begin()+i);
                }
            }
        }
    }

    for(const auto & v: reached_vertices){
        dij_distances[v] = null_distance;
        predecessors[v] = null_vertex;
        pruning_flag[v] = false;
    }
    reached_vertices.clear();
    L++;
}

void HighwayLabelling::AddLandmarkDirected(vertex r) {
    // compute Highway distances from previous landmarks to r
    // add distances in Highway between r and landmarks directly connected to r (i.e. in L(r))
    landmarks.insert(r);
    highway[r] = std::unordered_map<vertex, dist> ();
    for(size_t i = 0; i < out_landmarks_distances[r].size(); i++){
        highway[out_landmarks_distances[r][i]][r] = out_distances[r][i];
    }
    for(size_t i = 0; i < in_landmarks_distances[r].size(); i++){
        highway[r][in_landmarks_distances[r][i]] = in_distances[r][i];
    }

    highway[r][r] = 0;
    // fill distances in Highway between r and landmarks not covering r (i.e. not in L(r))
    for(const vertex & l: landmarks){
        if(highway[r].find(l) == highway[r].end()){
            highway[r][l] = null_distance;
            for(size_t j = 0; j < in_landmarks_distances[r].size(); j++){
                highway[r][l] = std::min(highway[r][l],highway[r][in_landmarks_distances[r][j]] + highway[in_landmarks_distances[r][j]][l]);
            }
        }
        if(highway[l].find(r) == highway[l].end()) {
            highway[l][r] = null_distance;
            for(size_t j = 0; j < out_landmarks_distances[r].size(); j++){
                highway[l][r] = std::min(highway[l][r],highway[out_landmarks_distances[r][j]][r] + highway[l][out_landmarks_distances[r][j]]);
            }
        }
    }

    // FORWARD

    out_landmarks_distances[r].clear();
    out_landmarks_distances[r].push_back(r);
    out_distances[r].clear();
    out_distances[r].push_back(0);
    // search from r to fill L a la Akiba (prune when reaching landmarks different from r, compute query to vertices not in L)
    std::vector<vertex> reached_vertices;
    std::vector<vertex> affected_vertices;
    dij_distances[r] = 0;
    reached_vertices.push_back(r);
    std::vector<vertex> reached_landmarks;
    std::priority_queue<std::pair<dist,vertex>, std::vector<std::pair<dist,vertex>>,
            PQComparator> pq;
    for(auto w: graph.neighborRange(r)){
        dij_distances[w] = dij_distances[r] + (dist)graph.weight(r,w);
        pq.push(std::make_pair(dij_distances[w], w));
        reached_vertices.push_back(w);
    }
    while(!pq.empty()){
        vertex v = pq.top().second;
        pq.pop();
        if(pruning_flag[v]) continue;
        if(landmarks.find(v) != landmarks.end() && v != r){
            pruning_flag[v] = true;
            reached_landmarks.push_back(v);
            continue;
        }
        dist query_dist = DirectedQueryDistance(r,v);
        if(query_dist < dij_distances[v]){
            pruning_flag[v] = true;
            continue;
        }
        affected_vertices.push_back(v);
        if(out_landmarks_distances[v].empty() || *out_landmarks_distances[v].rbegin() < r) {
            out_landmarks_distances[v].push_back(r);
            out_distances[v].push_back(dij_distances[v]);
        }
        else {
            for(size_t i = 0; i < out_landmarks_distances[v].size(); i++) {
                if(out_landmarks_distances[v][i] < r) {
                    continue;
                }
                if(out_landmarks_distances[v][i] > r) {
                    out_landmarks_distances[v].insert(out_landmarks_distances[v].begin()+i, r);
                    out_distances[v].insert(out_distances[v].begin()+i, dij_distances[v]);
                    break;
                }
                if(out_landmarks_distances[v][i] == r) {
                    out_distances[v][i] = dij_distances[v];
                    break;
                }
            }
        }

        for(auto w: graph.neighborRange(v)){
            if(dij_distances[w] > dij_distances[v] + graph.weight(v,w)){
                dij_distances[w] = dij_distances[v] + (dist)graph.weight(v,w);
                predecessors[w] = v;
                pruning_flag[w] = pruning_flag[v];
                pq.push(std::make_pair(dij_distances[w], w));
                reached_vertices.push_back(w);
            }

        }
    }

    for(const auto & v: reached_vertices){
        predecessors[v] = null_vertex;
        pruning_flag[v] = false;
    }


    // REVERSE
    in_landmarks_distances[r].clear();
    in_landmarks_distances[r].push_back(r);
    in_distances[r].clear();
    in_distances[r].push_back(0);
    // search from r to fill L a la Akiba (prune when reaching landmarks different from r, compute query to vertices not in L)
    dij_distances_in[r] = 0;
    std::vector<vertex> reached_vertices_in;
    std::vector<vertex> affected_vertices_in;
    std::vector<vertex> reached_landmarks_in;
    reached_vertices_in.push_back(r);
    while(!pq.empty()) pq.pop();
    for(auto w: graph.inNeighborRange(r)){
        dij_distances_in[w] = dij_distances_in[r] + (dist)graph.weight(w,r);
        pq.push(std::make_pair(dij_distances_in[w], w));
        reached_vertices_in.push_back(w);
    }
    while(!pq.empty()){
        vertex v = pq.top().second;
        pq.pop();
        if(pruning_flag[v]) continue;
        if(landmarks.find(v) != landmarks.end() && v != r){
            pruning_flag[v] = true;
            reached_landmarks_in.push_back(v);
            continue;
        }
        dist query_dist = DirectedQueryDistance(v, r);
        if(query_dist < dij_distances_in[v]){
            pruning_flag[v] = true;
            continue;
        }
        affected_vertices_in.push_back(v);
        if(in_landmarks_distances[v].empty() || *in_landmarks_distances[v].rbegin() < r) {
            in_landmarks_distances[v].push_back(r);
            in_distances[v].push_back(dij_distances_in[v]);
        }
        else {
            for(size_t i = 0; i < in_landmarks_distances[v].size(); i++) {
                if(in_landmarks_distances[v][i] < r) {
                    continue;
                }
                if(in_landmarks_distances[v][i] > r) {
                    in_landmarks_distances[v].insert(in_landmarks_distances[v].begin()+i, r);
                    in_distances[v].insert(in_distances[v].begin()+i, dij_distances_in[v]);
                    break;
                }
                if(in_landmarks_distances[v][i] == r) {
                    in_distances[v][i] = dij_distances_in[v];
                    break;
                }
            }
        }

        for(auto w: graph.inNeighborRange(v)){
            if(dij_distances_in[w] > dij_distances_in[v] + graph.weight(w, v)){
                dij_distances_in[w] = dij_distances_in[v] + (dist)graph.weight(w, v);
                predecessors[w] = v;
                pruning_flag[w] = pruning_flag[v];
                pq.push(std::make_pair(dij_distances_in[w], w));
                reached_vertices_in.push_back(w);
            }

        }
    }

    for(const auto & lndm: reached_landmarks){
        if(lndm == r) continue;
        for(const auto & v: affected_vertices_in){
            for(size_t i = 0; i < in_landmarks_distances[v].size(); i++){
                if(in_landmarks_distances[v][i] < lndm) continue;
                if(in_landmarks_distances[v][i] > lndm) break;
                if(in_distances[v][i] >= dij_distances_in[v] + highway[r][lndm]){
                    in_landmarks_distances[v].erase(in_landmarks_distances[v].begin()+i);
                    in_distances[v].erase(in_distances[v].begin()+i);
                }
            }
        }
    }

    // TODO: CHECK IF NEED AND CONTROL IN LINE 1241
    for(const auto & lndm: reached_landmarks_in){
        if(lndm == r) continue;
        for(const auto & v: affected_vertices){
            for(size_t i = 0; i < out_landmarks_distances[v].size(); i++){
                if(out_landmarks_distances[v][i] < lndm) continue;
                if(out_landmarks_distances[v][i] > lndm) break;
                if(out_distances[v][i] >= dij_distances[v] + highway[lndm][r]){
                    out_landmarks_distances[v].erase(out_landmarks_distances[v].begin()+i);
                    out_distances[v].erase(out_distances[v].begin()+i);
                }
            }
        }
    }

    for(const auto & v: reached_vertices){
        dij_distances[v] = null_distance;
    }

    for(const auto & v: reached_vertices_in){
        dij_distances_in[v] = null_distance;
        predecessors[v] = null_vertex;
        pruning_flag[v] = false;
    }
    reached_landmarks.clear();
    reached_landmarks_in.clear();
    reached_vertices.clear();
    reached_vertices_in.clear();
    L++;
}

void HighwayLabelling::RemoveLandmark(vertex r) {

    L--;
    landmarks.erase(r);

    std::vector<std::pair<vertex,dist>> affected_landmarks;

    // clear label of r and treat it as affected
    landmarks_distances[r].clear();
    distances[r].clear();
    // search from r for affected landmarks

    std::vector<vertex> reached_vertices;
    dij_distances[r] = 0;
    reached_vertices.push_back(r);
    std::priority_queue<std::pair<dist, vertex>, std::vector<std::pair<dist, vertex>>,
            PQComparator> pq;
    pq.push(std::make_pair(0, r));
    while (!pq.empty()) {
        vertex v = pq.top().second;
        pq.pop();
        if (pruning_flag[v]) continue;
        if (landmarks.find(v) != landmarks.end()) {
            pruning_flag[v] = true;
            if(dij_distances[v] > highway[r][v]) continue;
            affected_landmarks.emplace_back(v,dij_distances[v]);
            auto insertion_index = std::upper_bound(landmarks_distances[r].begin(),
                                                    landmarks_distances[r].end(), v);
            distances[r].insert(
                    distances[r].begin() + (insertion_index - landmarks_distances[r].begin()),
                    dij_distances[v]);
            landmarks_distances[r].insert(insertion_index, v);

            continue;
        }
        for(size_t i = 0; i < landmarks_distances[v].size(); i++){
            if(landmarks_distances[v][i] == r){
                landmarks_distances[v].erase(landmarks_distances[v].begin()+i);
                distances[v].erase(distances[v].begin()+i);
                break;
            }
            if(landmarks_distances[v][i] > r){
                break;
            }

        }
        for (auto w: graph.neighborRange(v)) {
            if (dij_distances[w] > dij_distances[v] + graph.weight(v, w)) {
                dij_distances[w] = dij_distances[v] + (dist) graph.weight(v, w);
                predecessors[w] = v;
                pruning_flag[w] = pruning_flag[v];
                pq.push(std::make_pair(dij_distances[w], w));
                reached_vertices.push_back(w);
            }

        }
    }

    for(const auto & v: reached_vertices){
        dij_distances[v] = null_distance;
        predecessors[v] = null_vertex;
        pruning_flag[v] = false;
    }
    reached_vertices.clear();

    // remove entries of r in highway
    for(const vertex & l: landmarks){
        highway[l].erase(r);
    }
    highway[r].clear();
    highway.erase(r);

    // search from r rooted in each affected landmark to potentially cover vertices
    for(const auto & ld: affected_landmarks){
        const vertex &l = ld.first;
        const dist &d = ld.second;
        dij_distances[l] = 0;
        dij_distances[r] = d;
        reached_vertices.push_back(r);
        std::priority_queue<std::pair<dist, vertex>, std::vector<std::pair<dist, vertex>>,
                PQComparator> pq;
        for (auto w: graph.neighborRange(r)) {
            if (w == l) continue;
            dij_distances[w] = dij_distances[r] + (dist) graph.weight(r, w);
            reached_vertices.push_back(w);
            pq.push(std::make_pair(dij_distances[w], w));
        }
        pruning_flag[l] = true;
        pruning_flag[r] = true;
        reached_vertices.push_back(l);

        while (!pq.empty()) {
            vertex v = pq.top().second;
            pq.pop();
            if (pruning_flag[v]) continue;
            if (landmarks.find(v) != landmarks.end()) {
                pruning_flag[v] = true;
                continue;
            }
            if(QueryDistance(v,l) <= dij_distances[v]){
                pruning_flag[v] = true;
                continue;
            }
            auto insertion_index = std::upper_bound(landmarks_distances[v].begin(),
                                               landmarks_distances[v].end(), l);
            distances[v].insert(
                    distances[v].begin() + (insertion_index - landmarks_distances[v].begin()),
                    dij_distances[v]);
            landmarks_distances[v].insert(insertion_index, l);


            for (auto w: graph.neighborRange(v)) {
                if (dij_distances[w] > dij_distances[v] + graph.weight(v, w)) {
                    dij_distances[w] = dij_distances[v] + (dist) graph.weight(v, w);
                    predecessors[w] = v;
                    pq.push(std::make_pair(dij_distances[w], w));
                    reached_vertices.push_back(w);
                }

            }
        }

        for(const auto & v: reached_vertices){
            dij_distances[v] = null_distance;
            predecessors[v] = null_vertex;
            pruning_flag[v] = false;
        }
        reached_vertices.clear();

    }
}


void HighwayLabelling::RemoveLandmarkDirected(vertex r) {

    L--;
    landmarks.erase(r);

    std::vector<std::pair<vertex,dist>> in_affected_landmarks;

    // clear label of r and treat it as affected
    in_landmarks_distances[r].clear();
    in_distances[r].clear();
    // search from r for affected landmarks

    std::vector<vertex> in_reached_vertices;
    dij_distances[r] = 0;
    in_reached_vertices.push_back(r);
    std::priority_queue<std::pair<dist, vertex>, std::vector<std::pair<dist, vertex>>,
            PQComparator> pq;
    pq.push(std::make_pair(0, r));
    while (!pq.empty()) {
        vertex v = pq.top().second;
        pq.pop();
        if (pruning_flag[v]) continue;
        if (landmarks.find(v) != landmarks.end()) {
            pruning_flag[v] = true;
            if(dij_distances[v] > highway[r][v]) continue;
            in_affected_landmarks.emplace_back(v,dij_distances[v]);
            auto insertion_index = std::upper_bound(in_landmarks_distances[r].begin(),
                                                    in_landmarks_distances[r].end(), v);
            in_distances[r].insert(
                    in_distances[r].begin() + (insertion_index - in_landmarks_distances[r].begin()),
                    dij_distances[v]);
            in_landmarks_distances[r].insert(insertion_index, v);

            continue;
        }
        for(size_t i = 0; i < out_landmarks_distances[v].size(); i++){
            if(out_landmarks_distances[v][i] == r){
                out_landmarks_distances[v].erase(out_landmarks_distances[v].begin()+i);
                out_distances[v].erase(out_distances[v].begin()+i);
                break;
            }
            if(out_landmarks_distances[v][i] > r){
                break;
            }

        }
        for (auto w: graph.neighborRange(v)) {
            if (dij_distances[w] > dij_distances[v] + graph.weight(v, w)) {
                dij_distances[w] = dij_distances[v] + (dist) graph.weight(v, w);
                predecessors[w] = v;
                pruning_flag[w] = pruning_flag[v];
                pq.push(std::make_pair(dij_distances[w], w));
                in_reached_vertices.push_back(w);
            }

        }
    }

    for(const auto & v: in_reached_vertices){
        dij_distances[v] = null_distance;
        predecessors[v] = null_vertex;
        pruning_flag[v] = false;
    }

    std::vector<std::pair<vertex,dist>> out_affected_landmarks;
    std::vector<vertex> out_reached_vertices;

    // clear label of r and treat it as affected
    out_landmarks_distances[r].clear();
    out_distances[r].clear();
    // search from r for affected landmarks

    dij_distances[r] = 0;
    out_reached_vertices.push_back(r);
    while(!pq.empty()) pq.pop();
    pq.push(std::make_pair(0, r));
    while (!pq.empty()) {
        vertex v = pq.top().second;
        pq.pop();
        if (pruning_flag[v]) continue;
        if (landmarks.find(v) != landmarks.end()) {
            pruning_flag[v] = true;
            if(dij_distances[v] > highway[v][r]) continue;
            out_affected_landmarks.emplace_back(v,dij_distances[v]);
            auto insertion_index = std::upper_bound(out_landmarks_distances[r].begin(),
                                                    out_landmarks_distances[r].end(), v);
            out_distances[r].insert(
                    out_distances[r].begin() + (insertion_index - out_landmarks_distances[r].begin()),
                    dij_distances[v]);
            out_landmarks_distances[r].insert(insertion_index, v);

            continue;
        }
        for(size_t i = 0; i < in_landmarks_distances[v].size(); i++){
            if(in_landmarks_distances[v][i] == r){
                in_landmarks_distances[v].erase(in_landmarks_distances[v].begin()+i);
                in_distances[v].erase(in_distances[v].begin()+i);
                break;
            }
            if(in_landmarks_distances[v][i] > r){
                break;
            }

        }
        for (auto w: graph.inNeighborRange(v)) {
            if (dij_distances[w] > dij_distances[v] + graph.weight(w, v)) {
                dij_distances[w] = dij_distances[v] + (dist) graph.weight(w, v);
                predecessors[w] = v;
                pruning_flag[w] = pruning_flag[v];
                pq.push(std::make_pair(dij_distances[w], w));
                out_reached_vertices.push_back(w);
            }

        }
    }

    for(const auto & v: out_reached_vertices){
        dij_distances[v] = null_distance;
        predecessors[v] = null_vertex;
        pruning_flag[v] = false;
    }

    // remove entries of r in highway
    for(const vertex & l: landmarks){
        highway[l].erase(r);
    }
    highway[r].clear();
    highway.erase(r);

    // search from r rooted in each affected landmark to potentially cover vertices
    for(const auto & ld: in_affected_landmarks){
        const vertex &l = ld.first;
        const dist &d = ld.second;
        dij_distances[l] = 0;
        dij_distances[r] = d;
        in_reached_vertices.push_back(r);
        std::priority_queue<std::pair<dist, vertex>, std::vector<std::pair<dist, vertex>>,
                PQComparator> pq;
        for (auto w: graph.inNeighborRange(r)) {
            if (w == l) continue;
            dij_distances[w] = dij_distances[r] + (dist) graph.weight(w, r);
            in_reached_vertices.push_back(w);
            pq.push(std::make_pair(dij_distances[w], w));
        }
        pruning_flag[l] = true;
        pruning_flag[r] = true;
        in_reached_vertices.push_back(l);

        while (!pq.empty()) {
            vertex v = pq.top().second;
            pq.pop();
            if (pruning_flag[v]) continue;
            if (landmarks.find(v) != landmarks.end()) {
                pruning_flag[v] = true;
                continue;
            }
            if(DirectedQueryDistance(v,l) < dij_distances[v]){
                pruning_flag[v] = true;
                continue;
            }
            if(*in_landmarks_distances[v].rbegin() < l || in_landmarks_distances.empty()) {
                in_landmarks_distances[v].push_back(l);
                in_distances[v].push_back(dij_distances[v]);
            }
            else {
                for(size_t i = 0; i < in_landmarks_distances[v].size(); i++) {
                    if(in_landmarks_distances[v][i] < l) {
                        continue;
                    }
                    if(in_landmarks_distances[v][i] > l) {
                        in_distances[v].insert(in_distances[v].begin() + i, dij_distances[v]);
                        in_landmarks_distances[v].insert(in_landmarks_distances[v].begin() + i, l);
                        break;
                    }
                    if(in_landmarks_distances[v][i] == l) {
                        in_distances[v][i] = dij_distances[v];
                        break;
                    }
                }
            }

            for (auto w: graph.inNeighborRange(v)) {
                if (dij_distances[w] > dij_distances[v] + graph.weight(w, v)) {
                    dij_distances[w] = dij_distances[v] + (dist) graph.weight(w, v);
                    predecessors[w] = v;
                    pq.push(std::make_pair(dij_distances[w], w));
                    in_reached_vertices.push_back(w);
                }

            }
        }

        for(const auto & v: in_reached_vertices){
            dij_distances[v] = null_distance;
            predecessors[v] = null_vertex;
            pruning_flag[v] = false;
        }
        in_reached_vertices.clear();

    }

    // search from r rooted in each affected landmark to potentially cover vertices
    for(const auto & ld: out_affected_landmarks){
        const vertex &l = ld.first;
        const dist &d = ld.second;
        dij_distances[l] = 0;
        dij_distances[r] = d;
        out_reached_vertices.push_back(r);
        std::priority_queue<std::pair<dist, vertex>, std::vector<std::pair<dist, vertex>>,
                PQComparator> pq;
        for (auto w: graph.neighborRange(r)) {
            if (w == l) continue;
            dij_distances[w] = dij_distances[r] + (dist) graph.weight(r, w);
            out_reached_vertices.push_back(w);
            pq.push(std::make_pair(dij_distances[w], w));
        }
        pruning_flag[l] = true;
        pruning_flag[r] = true;
        out_reached_vertices.push_back(l);

        while (!pq.empty()) {
            vertex v = pq.top().second;
            pq.pop();
            if (pruning_flag[v]) continue;
            if (landmarks.find(v) != landmarks.end()) {
                pruning_flag[v] = true;
                continue;
            }

            if(DirectedQueryDistance(l,v) < dij_distances[v]){
                pruning_flag[v] = true;
                continue;
            }

            if(out_landmarks_distances[v].empty() || *out_landmarks_distances[v].rbegin() < l) {
                    out_landmarks_distances[v].push_back(l);
                    out_distances[v].push_back(dij_distances[v]);
            }
            for(size_t i = 0; i < out_landmarks_distances[v].size(); i++) {
                if(out_landmarks_distances[v][i] < l) {
                    continue;
                }
                if(out_landmarks_distances[v][i] > l) {
                    out_distances[v].insert(out_distances[v].begin() + i, dij_distances[v]);
                    out_landmarks_distances[v].insert(out_landmarks_distances[v].begin() + i, l);
                    break;
                }
                if(out_landmarks_distances[v][i] == l) {
                    out_distances[v][i] = dij_distances[v];
                    break;
                }
            }

            for (auto w: graph.neighborRange(v)) {
                if (dij_distances[w] > dij_distances[v] + graph.weight(v, w)) {
                    dij_distances[w] = dij_distances[v] + (dist) graph.weight(v, w);
                    predecessors[w] = v;
                    pq.push(std::make_pair(dij_distances[w], w));
                    out_reached_vertices.push_back(w);
                }

            }
        }

        for(const auto & v: out_reached_vertices){
            dij_distances[v] = null_distance;
            predecessors[v] = null_vertex;
            pruning_flag[v] = false;
        }
        out_reached_vertices.clear();

    }

}

dist HighwayLabelling::BoundedSearch(vertex s, vertex t) {
    // Initialization
    dist *to_vertices = new dist[V];
    for(vertex j = 0; j < V; j++) {
        to_vertices[j] = null_distance;
    }
    dist m = null_distance;
    for(vertex v = 0; v < L; v++){
        dist s_to_v = P2PBFS(s, reverse_ordering[v], m, 0, to_vertices);
        if(s_to_v == null_distance)
            continue;
        dist s_to_t = P2PBFS(t, reverse_ordering[v], m, s_to_v, to_vertices);
        m = min(m, s_to_t);
    }

    delete [] to_vertices;
    return m;
}

inline void HighwayLabelling::StoreIndex(std::string filename, bool dynamic) {
        std::ofstream ofs(std::string("index/")+std::string(filename) +
            (dynamic ? std::string("_dynamic_") : std::string("_scratch_"))+ std::string("index"));
    for (int i = 0; i < V; i++) {
        vertex C = 0;
        for (int j = 0; j < distances[i].size(); j++) {
            if(distances[i][j] != null_distance)
                C++;
        }
        ofs.write((char*)&C, sizeof(C));
        for (int j = 0; j < distances[i].size(); j++) {
            if(distances[i][j] != null_distance) {
                ofs.write((char*)&j, sizeof(j));
                ofs.write((char*)&distances[i][j], sizeof(distances[i][j]));
            }
        }
    }

    for (const vertex & v: landmarks) {
        for (const vertex & w: landmarks) {
            if(highway[v][w] != null_distance)
                ofs.write((char*)&highway[v][w], sizeof(highway[v][w]));
        }
    }
    ofs.close();
}


#endif  // HGHWAY_LABELING_H_
