// ---- simple_prim.hpp ----
#pragma once
#include <vector>
#include <limits>
#include <cassert>

// lightweight undirected weighted graph built from edge list + weights
struct prim_graph {
    int n;
    std::vector<float> W;    

    // edges: array of (u,v) pairs; weights: array of float; m: #edges
    prim_graph(int n_, const std::pair<int,int>* edges, const float* weights, int m) : n(n_), W(size_t(n_) * size_t(n_), std::numeric_limits<float>::infinity()) {
        // zero diagonal
        for (int i = 0; i < n; ++i) W[size_t(i)*n + i] = 0.0f;
        
        for (int k = 0; k < m; ++k) {
            int u = edges[k].first;
            int v = edges[k].second;
            float w = weights[k];
            assert(0 <= u && u < n && 0 <= v && v < n);
            W[size_t(u)*n + v] = w;
            W[size_t(v)*n + u] = w;
        }
    }

    int num_vertices() const { return n; }
    
    float weight(int u, int v) const { return W[size_t(u)*n + v]; }
};

// Replacement for prim_minimum_spanning_tree(g, &p[0]) in boost
inline void prim(const prim_graph& g, int* parent_out) {
    const int n = g.num_vertices();
    assert(n > 0);

    
    const float INF = std::numeric_limits<float>::infinity();
    std::vector<float> key(n, INF);
    std::vector<char>  inMST(n, 0);

    // root = 0
    key[0] = 0.0f;
    parent_out[0] = 0;

    for (int it = 0; it < n; ++it) {
        int u = -1;
        float best = INF;
        for (int v = 0; v < n; ++v) {
            if (!inMST[v] && key[v] < best) {
                best = key[v];
                u = v;
            }
        }
        if (u == -1) break; 
        inMST[u] = 1;
        
        for (int v = 0; v < n; ++v) {
            if (inMST[v] || v == u) continue;
            float w = g.weight(u, v);
            if (w < key[v]) {
                key[v] = w;
                parent_out[v] = u;
            }
        }
    }
}
