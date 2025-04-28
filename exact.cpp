#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <limits>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>

using namespace std;

class Graph {
private:
    int n;
    vector<vector<int>> adj;
    
public:
    Graph(int vertices) : n(vertices) {
        adj.resize(n);
    }
    
    Graph() : n(0) {}
    
    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    const vector<vector<int>>& getAdjList() const {
        return adj;
    }
    
    int getVertexCount() const {
        return n;
    }
    
    vector<int> computeCoreness() const {
        vector<int> degree(n);
        vector<int> core(n);
        vector<vector<int>> bin(n + 1);
        
        for (int v = 0; v < n; v++) {
            degree[v] = adj[v].size();
            bin[degree[v]].push_back(v);
        }
        
        int maxDegree = 0;
        for (int v = 0; v < n; v++) {
            maxDegree = max(maxDegree, degree[v]);
        }
        
        vector<bool> processed(n, false);
        for (int d = 0; d <= maxDegree; d++) {
            while (!bin[d].empty()) {
                int v = bin[d].back();
                bin[d].pop_back();
                processed[v] = true;
                core[v] = d;
                
                for (int u : adj[v]) {
                    if (!processed[u] && degree[u] > d) {
                        bin[degree[u]].erase(
                            find(bin[degree[u]].begin(), bin[degree[u]].end(), u)
                        );
                        degree[u]--;
                        bin[degree[u]].push_back(u);
                    }
                }
            }
        }
        
        return core;
    }
    
    vector<vector<int>> findTriangles() const {
        vector<vector<int>> triangles;
        
        for (int i = 0; i < n; i++) {
            unordered_set<int> neighbors(adj[i].begin(), adj[i].end());
            
            for (int j : adj[i]) {
                if (j > i) {
                    for (int k : adj[j]) {
                        if (k > j && neighbors.count(k)) {
                            triangles.push_back({i, j, k});
                        }
                    }
                }
            }
        }
        
        return triangles;
    }
    
    int countTriangles() const {
        int count = 0;
        for (int i = 0; i < n; i++) {
            unordered_set<int> neighbors(adj[i].begin(), adj[i].end());
            
            for (int j : adj[i]) {
                if (j > i) {
                    for (int k : adj[j]) {
                        if (k > j && neighbors.count(k)) {
                            count++;
                        }
                    }
                }
            }
        }
        return count;
    }
    
    vector<vector<int>> findHMinus1Cliques(int h) const {
        if (h == 2) {
            vector<vector<int>> edges;
            for (int i = 0; i < n; i++) {
                for (int j : adj[i]) {
                    if (i < j) { 
                        edges.push_back({i, j});
                    }
                }
            }
            return edges;
        }
        
        if (h == 3) {
            vector<vector<int>> edges;
            for (int i = 0; i < n; i++) {
                for (int j : adj[i]) {
                    if (i < j) { 
                        edges.push_back({i, j});
                    }
                }
            }
            return edges;
        }
        
        vector<vector<int>> cliques;
        vector<int> currentClique;
        
        function<void(int, int)> findCliques = [&](int start, int depth) {
            if (depth == h - 1) {
                cliques.push_back(currentClique);
                return;
            }
            
            for (int i = start; i < n; i++) {
                bool canAdd = true;
                
                for (int v : currentClique) {
                    bool isConnected = false;
                    for (int neighbor : adj[v]) {
                        if (neighbor == i) {
                            isConnected = true;
                            break;
                        }
                    }
                    if (!isConnected) {
                        canAdd = false;
                        break;
                    }
                }
                
                if (canAdd) {
                    currentClique.push_back(i);
                    findCliques(i + 1, depth + 1);
                    currentClique.pop_back();
                }
            }
        };
        
        findCliques(0, 0);
        return cliques;
    }

    void printSubgraphEdges(ofstream& outFile) const {
        for (int i = 0; i < n; i++) {
            for (int j : adj[i]) {
                if (i < j) {
                    outFile << "(" << (i+1) << ", " << (j+1) << ")" << endl;
                }
            }
        }
    }    
    
    int cliqueDegreee(int v, int h) const {
        if (h == 2) {
            return adj[v].size();
        }
        
        if (h == 3) {
            int count = 0;
            unordered_set<int> neighbors(adj[v].begin(), adj[v].end());
            
            for (int u : adj[v]) {
                for (int w : adj[u]) {
                    if (w != v && neighbors.count(w)) {
                        count++;
                    }
                }
            }
            return count / 2;
        }
        
        int degree = 0;
        
        vector<int> vertices;
        for (int i = 0; i < n; i++) {
            if (i != v) vertices.push_back(i);
        }
        
        function<void(int, int, vector<int>&)> countCliques = [&](int start, int depth, vector<int>& current) {
            if (depth == h - 1) {
                bool isClique = true;
                for (int u : current) {
                    bool connected = false;
                    for (int neighbor : adj[v]) {
                        if (neighbor == u) {
                            connected = true;
                            break;
                        }
                    }
                    if (!connected) {
                        isClique = false;
                        break;
                    }
                }
                
                for (int i = 0; i < current.size(); i++) {
                    for (int j = i + 1; j < current.size(); j++) {
                        bool connected = false;
                        for (int neighbor : adj[current[i]]) {
                            if (neighbor == current[j]) {
                                connected = true;
                                break;
                            }
                        }
                        if (!connected) {
                            isClique = false;
                            break;
                        }
                    }
                    if (!isClique) break;
                }
                
                if (isClique) degree++;
                return;
            }
            
            for (int i = start; i < vertices.size(); i++) {
                current.push_back(vertices[i]);
                countCliques(i + 1, depth + 1, current);
                current.pop_back();
            }
        };
        
        vector<int> current;
        countCliques(0, 0, current);
        
        return degree;
    }
    
    bool formsHClique(int v, const vector<int>& clique) const {
        for (int u : clique) {
            bool connected = false;
            for (int neighbor : adj[v]) {
                if (neighbor == u) {
                    connected = true;
                    break;
                }
            }
            if (!connected) return false;
        }
        return true;
    }
    
    Graph extractSubgraph(const vector<int>& vertices) const {
        unordered_map<int, int> vertexMap;
        for (int i = 0; i < vertices.size(); i++) {
            vertexMap[vertices[i]] = i;
        }
    
        Graph subgraph(vertices.size());
        for (int i = 0; i < vertices.size(); i++) {
            int v = vertices[i];
            for (int u : adj[v]) {
                auto it = vertexMap.find(u);
                if (it != vertexMap.end() && i < it->second) {
                    subgraph.addEdge(i, it->second);
                }
            }
        }
        return subgraph;
    }
    
    double calculateHCliqueDensity(int h) const {
        if (n <= 1) return 0.0;
        
        int numCliques = 0;
        
        if (h == 2) {
            for (int i = 0; i < n; i++) {
                for (int j : adj[i]) {
                    if (i < j) {
                        numCliques++;
                    }
                }
            }
        } else if (h == 3) {
            numCliques = countTriangles();
        } else {
            vector<vector<int>> hMinus1Cliques = findHMinus1Cliques(h);
            for (int v = 0; v < n; v++) {
                for (const auto& clique : hMinus1Cliques) {
                    if (formsHClique(v, clique)) {
                        numCliques++;
                    }
                }
            }
            numCliques /= h;
        }
        
        return (double)numCliques / n;
    }
};

class SparseMaxFlow {
private:
    int n;
    unordered_map<int, unordered_map<int, int>> capacity;
    unordered_map<int, unordered_map<int, int>> flow;
    
public:
    SparseMaxFlow(int vertices) : n(vertices) {}
    
    void addEdge(int u, int v, int cap) {
        capacity[u][v] = cap;
    }
    
    int fordFulkerson(int source, int sink) {
        for (auto& [u, edges] : capacity) {
            for (auto& [v, cap] : edges) {
                flow[u][v] = 0;
            }
        }
        
        int maxFlow = 0;
        
        while (true) {
            vector<int> parent(n, -1);
            queue<int> q;
            q.push(source);
            parent[source] = -2;
            
            while (!q.empty() && parent[sink] == -1) {
                int u = q.front();
                q.pop();
                
                for (auto& [v, cap] : capacity[u]) {
                    if (parent[v] == -1 && cap > flow[u][v]) {
                        parent[v] = u;
                        q.push(v);
                    }
                }
            }
            
            if (parent[sink] == -1) break;
            
            int pathFlow = numeric_limits<int>::max();
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                pathFlow = min(pathFlow, capacity[u][v] - flow[u][v]);
            }
            
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                flow[u][v] += pathFlow;
                flow[v][u] -= pathFlow;
            }
            
            maxFlow += pathFlow;
        }
        
        return maxFlow;
    }
    
    pair<vector<int>, vector<int>> minCut(int source, int sink) {
        fordFulkerson(source, sink);
        
        vector<bool> visited(n, false);
        queue<int> q;
        q.push(source);
        visited[source] = true;
        
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            
            for (auto& [v, cap] : capacity[u]) {
                if (!visited[v] && cap > flow[u][v]) {
                    visited[v] = true;
                    q.push(v);
                }
            }
        }
        
        vector<int> S, T;
        for (int i = 0; i < n; i++) {
            if (visited[i]) S.push_back(i);
            else T.push_back(i);
        }
        
        return {S, T};
    }
};

Graph readGraphFromFile(const string& filename) {
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        exit(1);
    }
    
    int n, e;
    file >> n >> e;
    Graph G(n);
    
    string line;
    getline(file, line);
    
    for (int i = 0; i < e; i++) {
        int from, to;
        file >> from >> to;
        from--;
        to--;
        G.addEdge(from, to);
    }
    
    file.close();
    return G;
}

Graph exactDensestSubgraph(const Graph& G, int h, ofstream& outFile);

Graph exactDensestSubgraphWithPruning(const Graph& G, int h, ofstream& outFile) {
    int n = G.getVertexCount();
    
    vector<int> cores = G.computeCoreness();
    int kmax = *max_element(cores.begin(), cores.end());
    
    double L = kmax / 2.0;
    
    vector<int> coreVertices;
    for (int v = 0; v < n; v++) {
        if (cores[v] >= ceil(L)) {
            coreVertices.push_back(v);
        }
    }
    
    if (coreVertices.empty() || coreVertices.size() == n) {
        return exactDensestSubgraph(G, h, outFile);
    }
    
    Graph coreGraph = G.extractSubgraph(coreVertices);
    
    Graph densestSubgraph = exactDensestSubgraph(coreGraph, h, outFile);
    
    vector<int> originalVertices;
    for (int i = 0; i < densestSubgraph.getVertexCount(); i++) {
        originalVertices.push_back(coreVertices[i]);
    }
    
    return G.extractSubgraph(originalVertices);
}

Graph exactDensestSubgraph(const Graph& G, int h, ofstream& outFile) {
    int n = G.getVertexCount();
    
    double l = 0;
    double u = 0;
    
    for (int v = 0; v < n; v++) {
        u = max(u, (double)G.cliqueDegreee(v, h));
    }
    
    vector<vector<int>> hMinus1Cliques = G.findHMinus1Cliques(h);
    outFile << "Found " << hMinus1Cliques.size() << " (" << (h-1) << ")-cliques" << endl;
    
    if (h == 3) {
        int triangleCount = G.countTriangles();
        outFile << "Found " << triangleCount << " triangles" << endl;
    }
    
    if (hMinus1Cliques.empty()) {
        return Graph(0);
    }
    
    vector<int> densestSubgraphVertices;
    
    while (u - l > 1e-6) {
        double alpha = (l + u) / 2;
        
        int networkSize = 1 + n + hMinus1Cliques.size() + 1;
        SparseMaxFlow flowNetwork(networkSize);
        
        int source = 0;
        int sink = networkSize - 1;
        int vertexOffset = 1;
        int cliqueOffset = vertexOffset + n;
        
        for (int v = 0; v < n; v++) {
            flowNetwork.addEdge(source, vertexOffset + v, G.cliqueDegreee(v, h));
        }
        
        for (int v = 0; v < n; v++) {
            flowNetwork.addEdge(vertexOffset + v, sink, alpha * h);
        }
        
        for (int i = 0; i < hMinus1Cliques.size(); i++) {
            for (int v : hMinus1Cliques[i]) {
                flowNetwork.addEdge(cliqueOffset + i, vertexOffset + v, numeric_limits<int>::max());
            }
        }
        
        for (int v = 0; v < n; v++) {
            for (int i = 0; i < hMinus1Cliques.size(); i++) {
                if (G.formsHClique(v, hMinus1Cliques[i])) {
                    flowNetwork.addEdge(vertexOffset + v, cliqueOffset + i, 1);
                }
            }
        }
        
        auto [S, T] = flowNetwork.minCut(source, sink);
        
        if (S.size() == 1 && S[0] == source) {
            u = alpha;
        } else {
            l = alpha;
            
            vector<int> subgraphVertices;
            for (int node : S) {
                if (node != source && node >= vertexOffset && node < cliqueOffset) {
                    subgraphVertices.push_back(node - vertexOffset);
                }
            }
            
            densestSubgraphVertices = subgraphVertices;
        }
    }
    
    if (densestSubgraphVertices.empty()) {
        return G;
    }
    
    return G.extractSubgraph(densestSubgraphVertices);
}

int main() {    
    auto start = chrono::high_resolution_clock::now();
    Graph G = readGraphFromFile("CA-HepTh-preprocessed.txt");
    
    int h;
    cout << "Enter H value: ";
    cin >> h;
    
    ofstream outFile("CA-HepTh-output.txt");
    
    Graph densestSubgraph;
    densestSubgraph = exactDensestSubgraph(G, h, outFile);
    
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    
    outFile << "Densest subgraph found with " << densestSubgraph.getVertexCount() << " vertices. Edges:\n\n";
    densestSubgraph.printSubgraphEdges(outFile);
    outFile << "\nComputation time: " << elapsed.count() << " seconds" << endl;
    
    outFile.close();
    
    return 0;
}