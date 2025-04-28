#include <bits/stdc++.h>
using namespace std;

const int MAXN = 2000;

struct Edge {
    int to;
    double cap, flow;
    int rev;
};

struct Graph {
    int n;
    vector<vector<int>> adj;
    Graph(int _n) : n(_n), adj(_n) {}
    void add_edge(int u, int v) {
        if (u == v || find(adj[u].begin(), adj[u].end(), v) != adj[u].end()) return;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    bool is_connected(int u, int v) const {
        return find(adj[u].begin(), adj[u].end(), v) != adj[u].end();
    }
};

struct FlowGraph {
    int n;
    vector<vector<Edge>> adj;
    FlowGraph(int _n) : n(_n), adj(_n) {}
    void add_edge(int from, int to, double cap) {
        adj[from].push_back({to, cap, 0, (int)adj[to].size()});
        adj[to].push_back({from, 0, 0, (int)adj[from].size() - 1});
    }
    double max_flow(int s, int t) {
        vector<int> level(n, -1);
        auto bfs = [&]() {
            fill(level.begin(), level.end(), -1);
            queue<int> q;
            level[s] = 0;
            q.push(s);
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                for (auto& e : adj[u]) {
                    if (level[e.to] == -1 && e.cap - e.flow > 1e-6) {
                        level[e.to] = level[u] + 1;
                        q.push(e.to);
                    }
                }
            }
            return level[t] != -1;
        };
        function<double(int, int, double)> dfs = [&](int u, int t, double f) {
            if (u == t) return f;
            double pushed = 0;
            for (auto& e : adj[u]) {
                if (level[e.to] == level[u] + 1 && level[e.to] > 0) {
                    if (e.cap - e.flow > 1e-6) {
                        double flow = min(f, e.cap - e.flow);
                        flow = dfs(e.to, t, flow);
                        if (flow > 1e-6) {
                            e.flow += flow;
                            adj[e.to][e.rev].flow -= flow;
                            pushed += flow;
                        }
                    }
                }
            }
            return pushed;
        };
        double flow = 0;
        while (bfs()) {
            flow += dfs(s, t, 1e9);
        }
        return flow;
    }
};

vector<int> find_min_cut(const FlowGraph& fg, int s) {
    vector<bool> visited(fg.n, false);
    queue<int> q;
    q.push(s);
    visited[s] = true;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (auto& e : fg.adj[u]) {
            if (!visited[e.to] && e.cap - e.flow > 1e-6) {
                visited[e.to] = true;
                q.push(e.to);
            }
        }
    }
    vector<int> S;
    for (int i = 0; i < fg.n; i++) {
        if (visited[i]) S.push_back(i);
    }
    return S;
}

vector<int> core_decomposition(const Graph& G) {
    vector<int> deg(G.n, 0);
    int max_deg = 0;
    for (int i = 0; i < G.n; i++) {
        deg[i] = G.adj[i].size();
        if (deg[i] > max_deg) max_deg = deg[i];
    }
    vector<vector<int>> bucket(max_deg + 1);
    for (int i = 0; i < G.n; i++) {
        bucket[deg[i]].push_back(i);
    }
    vector<bool> processed(G.n, false);
    vector<int> core(G.n);
    int k = 0;
    for (int d = 0; d <= max_deg; d++) {
        while (!bucket[d].empty()) {
            int v = bucket[d].back();
            bucket[d].pop_back();
            if (processed[v]) continue;
            processed[v] = true;
            core[v] = k;
            for (int u : G.adj[v]) {
                if (!processed[u]) {
                    int current_deg = deg[u];
                    deg[u]--;
                    if (deg[u] < current_deg) {
                        bucket[deg[u]].push_back(u);
                    }
                }
            }
        }
        k++;
    }
    return core;
}

Graph extract_core(const Graph& G, const vector<int>& core, int k) {
    vector<int> vertices;
    for (int i = 0; i < G.n; i++) {
        if (core[i] >= k) vertices.push_back(i);
    }
    Graph core_graph(vertices.size());
    map<int, int> old_to_new;
    for (int i = 0; i < vertices.size(); i++) {
        old_to_new[vertices[i]] = i;
    }
    for (int v : vertices) {
        for (int u : G.adj[v]) {
            if (core[u] >= k && u > v) {
                core_graph.add_edge(old_to_new[v], old_to_new[u]);
            }
        }
    }
    return core_graph;
}

vector<vector<int>> find_connected_components(const Graph& G) {
    vector<bool> visited(G.n, false);
    vector<vector<int>> components;
    for (int v = 0; v < G.n; v++) {
        if (!visited[v]) {
            vector<int> component;
            queue<int> q;
            q.push(v);
            visited[v] = true;
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                component.push_back(u);
                for (int w : G.adj[u]) {
                    if (!visited[w]) {
                        visited[w] = true;
                        q.push(w);
                    }
                }
            }
            components.push_back(component);
        }
    }
    return components;
}

vector<vector<int>> enumerate_triangles(const Graph& G) {
    vector<vector<int>> triangles;
    set<tuple<int, int, int>> tri_set;
    for (int v = 0; v < G.n; v++) {
        for (size_t i = 0; i < G.adj[v].size(); i++) {
            for (size_t j = i + 1; j < G.adj[v].size(); j++) {
                int u = G.adj[v][i];
                int w = G.adj[v][j];
                if (u < w && G.is_connected(u, w)) {
                    vector<int> tri = {u, v, w};
                    sort(tri.begin(), tri.end());
                    tri_set.insert({tri[0], tri[1], tri[2]});
                }
            }
        }
    }
    for (auto& t : tri_set) {
        triangles.push_back({get<0>(t), get<1>(t), get<2>(t)});
    }
    return triangles;
}

double compute_density(const Graph& G, const vector<int>& vertices) {
    if (vertices.empty()) return 0;
    set<int> vertex_set(vertices.begin(), vertices.end());
    int triangle_count = 0;
    for (int v : vertices) {
        for (size_t i = 0; i < G.adj[v].size(); i++) {
            for (size_t j = i + 1; j < G.adj[v].size(); j++) {
                int u = G.adj[v][i];
                int w = G.adj[v][j];
                if (vertex_set.count(u) && vertex_set.count(w) && G.is_connected(u, w)) {
                    triangle_count++;
                }
            }
        }
    }
    return triangle_count / (double)vertices.size();
}

vector<int> core_exact(const Graph& G, int h = 3) {
    vector<int> core = core_decomposition(G);
    int k_max = *max_element(core.begin(), core.end());
    double l = 0, u = k_max;
    vector<int> best_subgraph;
    double best_density = 0;
    const double epsilon = 1e-6;

    while (u - l > epsilon) {
        double g = (l + u) / 2;
        int k = ceil(g * h);
        if (k > k_max) k = k_max;
        Graph core_graph = extract_core(G, core, k);
        vector<vector<int>> components = find_connected_components(core_graph);
        bool feasible = false;

        for (const auto& comp_vertices : components) {
            if (comp_vertices.empty()) continue;
            Graph sub(comp_vertices.size());
            map<int, int> old_to_new;
            set<int> comp_set(comp_vertices.begin(), comp_vertices.end());
            for (int i = 0; i < comp_vertices.size(); i++) {
                old_to_new[comp_vertices[i]] = i;
            }
            for (int v : comp_vertices) {
                for (int u : G.adj[v]) {
                    if (core[u] >= k && comp_set.count(u) && u > v) {
                        sub.add_edge(old_to_new[v], old_to_new[u]);
                    }
                }
            }

            vector<vector<int>> triangles = enumerate_triangles(sub);
            int c = triangles.size();
            if (c == 0) continue;

            FlowGraph fg(2 + sub.n + c);
            for (int i = 0; i < c; i++) {
                fg.add_edge(0, sub.n + 2 + i, 1);
                for (int v : triangles[i]) {
                    fg.add_edge(sub.n + 2 + i, 2 + v, 1e9);
                }
            }
            for (int j = 0; j < sub.n; j++) {
                fg.add_edge(2 + j, 1, g);
            }

            double flow = fg.max_flow(0, 1);
            if (abs(flow - c) < 1e-6) {
                feasible = true;
            } else {
                vector<int> S = find_min_cut(fg, 0);
                vector<int> subgraph_vertices;
                for (int v : S) {
                    if (v >= 2 && v < 2 + sub.n) {
                        subgraph_vertices.push_back(comp_vertices[v - 2]);
                    }
                }
                double density = compute_density(G, subgraph_vertices);
                if (density > best_density) {
                    best_density = density;
                    best_subgraph = subgraph_vertices;
                }
            }
        }

        if (feasible) {
            l = g;
        } else {
            u = g;
        }
    }

    return best_subgraph;
}

int main() {
    ifstream infile("as19991212-preprocessed.txt");
    if (!infile.is_open()) {
        cerr << "Error: Unable to open dataset" << endl;
        return 1;
    }

    int n, m;
    infile >> n >> m;
    Graph G(n);
    for (int i = 0; i < m; i++) {
        int u, v;
        infile >> u >> v;
        u--;
        v--;
        if (u >= 0 && u < n && v >= 0 && v < n && u != v) {
            G.add_edge(u, v);
        }
    }
    infile.close();

    auto start = chrono::high_resolution_clock::now();

    int h = 3;
    cout << "Found " << m << " (" << (h-1) << ")-cliques" << endl;
    vector<vector<int>> triangles = enumerate_triangles(G);
    cout << "Found " << triangles.size() << " triangles" << endl;

    vector<int> densest_subgraph = core_exact(G);
    cout << "Densest subgraph found with " << densest_subgraph.size() << " vertices. Edges:" << endl;
    for (size_t i = 0; i < densest_subgraph.size(); i++) {
        for (size_t j = i + 1; j < densest_subgraph.size(); j++) {
            int u = densest_subgraph[i];
            int v = densest_subgraph[j];
        }
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "Time taken: " << duration.count() << " seconds" << endl;

    return 0;
}
