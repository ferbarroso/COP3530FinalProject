//
//  main.cpp
//  Project2
//
//  Created by Fernando Barroso on 4/24/16.
//  Copyright © 2016 ferbarroso. All rights reserved.
//

#include <iostream>
#include <vector>
#include <queue>

using std::cin;
using std::cout;
using std::string;
using std::min;
using std::endl;
using std::vector;
using std::priority_queue;
using std::queue;

// Initialize Classes
class Node;
class Edge;
class LinkedList;
class Graph;
class Dijkstra;

// Initialize Methods
int getMagiSequenceInt(int array[], int arraySize);
vector<int> getMagiSequenceVector(int array[], int arraySize);
int getEditDistance(string word1, int l1, string word2, int l2);

// This is a basic edge class that connects two nodes
class Edge {
    Node *origin;
    Node *destination;
    int weight;
    
public:
    Edge(Node *beginNode, Node *endNode, int w) {
        origin = beginNode;
        destination = endNode;
        weight = w;
    }
    
    Node* getOrigin() {
        return origin;
    }
    
    Node* getDestination() {
        return destination;
    }
    
    int getWeight() {
        return weight;
    }
};

// This is a basic node class
class Node {
    
    // list of every edge from this node
    vector<Edge*> edges;
    
public:
    string charm;
    int id;
    Node *next;
    Node *previous = NULL;
    vector<int> lis;
    double minDistance = std::numeric_limits<double>::infinity();
    vector<Edge*> pubEdges;
    
    Node(string x, int identity) {
        charm = x;
        id = identity;
        next = NULL;
    }
    Node(string x, int identity, Node *n) {
        charm = x;
        id = identity;
        next = n;
    }
    
    // Adds an edge object between this node and another
    void addEdge(Node *node, int weight) {
        Edge *myEdge = new Edge(this, node, weight);
        edges.push_back(myEdge);
        pubEdges = edges;
    }
    
    // This prints all the edges and their end destinations
    void printEdges() {
        cout << charm << ":" << endl;
        for (int i = 0; i < edges.size(); i++) {
            Edge *edge = edges[i];
            cout << "Dest: " << edge->getDestination()->charm << " with weight: " << edge->getWeight() << endl;
        }
        cout << endl;
    }
};

// This is a basic linked list class
class LinkedList {
    Node *head;
    Node *tail;
public:
    
    Node *publicHead;
    
    LinkedList() {
        head = NULL;
        publicHead = NULL;
        tail = NULL;
    }
    
    // Adds a node to the end of the linked list
    void appendNodeToTail(Node *node) {
        Node *p;
        
        if(head == NULL) {
            head = node;
            publicHead = head;
            tail = head;
        }
        else {
            p = tail;
            p->next = node;
            tail = p->next;
        }
    }
    
    // Prints the charm of each node in the linkedlist
    void printListOfValues() {
        Node *p = head;
        while(p != NULL) {
            cout << p->charm << " ";
            p = p->next;
        }
    }
    
};

// This is a basic graph class
class Graph {
    int numOfNodes;
    LinkedList *adjNodeList;
    
public:
    
    // Looks weird but is a dynamic array of node pointers
    // It contains all nodes in this particular graph
    Node* *nodes;
    
    // A public reference to the number of nodes so others can't mess with it
    int pubNumOfNodes;
    
    Graph(int numOfVertices) {
        this->numOfNodes = numOfVertices;
        pubNumOfNodes = numOfVertices;
        adjNodeList = new LinkedList[numOfVertices];
        nodes = new Node*[numOfVertices];
    }
    
    void addNodeToArray(string charm, int id, vector<int> lis) {
        Node *node = new Node(charm, id);
        node->lis = lis;
        nodes[id] = node;
    }
    
    void addEdgeUnidirectional(Node *v1, Node *v2, int weight) {
        v1->addEdge(v2, weight);
        adjNodeList[v1->id].appendNodeToTail(v2);
    }
    
    void addEdgeBidirectional(Node *v1, Node *v2, int weight) {
        v1->addEdge(v2, weight);
        adjNodeList[v1->id].appendNodeToTail(v2);
        
        v2->addEdge(v1, weight);
        adjNodeList[v2->id].appendNodeToTail(v1);
    }
    
    // This adds edges between all nodes
    void setupGraph() {
        for(int i = 0; i < numOfNodes; i++) {
            for(int j = i + 1; j < numOfNodes; j++) {
                addEdgeBidirectional(nodes[i], nodes[j], getEditDistance(nodes[i]->charm, (int)nodes[i]->charm.length(), nodes[j]->charm, (int)nodes[j]->charm.length()));
            }
        }
    }
    
    // This will return the node id that has a desired charm
    int findIDWithCharm(string charm) {
        Node *p = nodes[0];
        
        while (p != NULL && p->charm != charm) {
            p = p->next;
        }
        
        if (p != NULL) {
            return p->id;
        }
        else {
            return -1;
        }
    }
    
    // Prints all nodes that can be reached from all nodes
    void printAllNodes() {
        for (int i = 0; i < numOfNodes; i++) {
            cout << "From node " << i << ":\n";
            adjNodeList[i].printListOfValues();
            cout << "\n";
        }
    }
    
    // Prints all nodes that can be reached from one node
    void printNodeFrom(int id) {
        cout << "From node " << id << ":\n";
        adjNodeList[id].printListOfValues();
    }
    
    // This prints the nodes a node can see and the weight of each edge between them
    void printEdges() {
        for (int i = 0; i < numOfNodes; i++) {
            cout << "From node " << i << ":\n";
            nodes[i]->printEdges();
            cout << "\n";
        }
    }
    
    // This prints the lis of each node
    void printNodeLis(int id) {
        vector<int> lis = nodes[id]->lis;
        
        for (int i = 0; i < lis.size(); i++) {
            cout << lis[i] << " ";
        }
        cout << "\n";
    }
};

// A class for holding Dijkstra's algorithm
class Dijkstra {
    
public:
    
    // We need a reference to the same graph
    Graph graph = NULL;
    
    Dijkstra(Graph g) {
        graph = g;
    }
    
    // This does the Dijkstra algorithm (a modified BFS) to find the minimum cost path
    void computePaths(Node *source) {
        
        // We want a priority queue to have the minimum on top
        source->minDistance = 0;
        priority_queue<Node*> nodeQueue;
        
        // We want to push the starting node on
        nodeQueue.push(source);
        
        // We go until we reach every node in the single connected graph
        while (!nodeQueue.empty()) {
            Node *u = nodeQueue.top();
            nodeQueue.pop();
            
            // Visit each edge that u can reach
            for (Edge *e : u->pubEdges) {
                Node *v = e->getDestination();
                double weight = e->getWeight();
                double distanceThroughU = u->minDistance + weight;
                
                // If we get a lower distance and also have enough to do the incantation
                if (distanceThroughU < v->minDistance && u->lis.size() >= e->getWeight()) {
                    
                    priority_queue<Node*> tempQueue;
                    bool done = false;
                    
                    // We remove the next from our queue
                    while (!nodeQueue.empty()) {
                        if (nodeQueue.top() == v && !done) {
                            nodeQueue.pop();
                            done = true;
                        }
                        else {
                            Node *p = nodeQueue.top();
                            tempQueue.push(p);
                            nodeQueue.pop();
                        }
                    }
                    
                    // We make put the modified queue back to our nodequeue
                    nodeQueue = tempQueue;
                    
                    v->minDistance = distanceThroughU;
                    v->previous = u; //predecesor
                    nodeQueue.push(v);
                }
            }
        }
    }
    
    // This starts at our end node and follows its path back to the beginning
    vector<Node*> getShortestPathTo(Node *target) {
        vector<Node*> path;
        for (Node *vertex = target; vertex != NULL; vertex = vertex->previous) {
            path.push_back(vertex);
        }
        
        // We flip it around to get the path from the beginning and not from the end
        std::reverse(path.begin(), path.end());
        
        return path;
    }
    
    // This prints out the appropriate stuff
    void printPath(vector<Node*> vec, int start, int end) {
        
        if (graph.nodes[start] != graph.nodes[end] && vec.size() == 1 && vec[vec.size()-1] == graph.nodes[end] && vec[vec.size()-1]->previous == NULL) {
            cout << "IMPOSSIBLE" << endl;
        }
        else {
            int edgeSum = 0;
            int gemSum = 0;
            
            for (int i = 0; i < vec.size() - 1; i++) {
                
                // We get the edit distance between the two nodes
                int editD = getEditDistance(vec[i]->charm, (int)vec[i]->charm.size(), vec[i + 1]->charm, (int)vec[i + 1]->charm.size());
                edgeSum = edgeSum + editD;
                
                // We sum only the editD number from the longest increasing sequence
                for (int j = 0; j < editD; j++) {
                    gemSum = gemSum + vec[i]->lis[j];
                }
            }
            
            cout << edgeSum << " " << gemSum << endl;
        }
    }
    
    // We need to reset the nodes to having infinity minDistance and no previous node pointer
    void resetNodes() {
        for (int i = 0; i < graph.pubNumOfNodes; i++) {
            graph.nodes[i]->previous = NULL;
            graph.nodes[i]->minDistance = std::numeric_limits<double>::infinity();
        }
    }
    
    // This is our controller method we call to do everything based on the paths
    void doDijkstra(int start, int end) {
        if (start == -1 || end == -1) {
            cout << "IMPOSSIBLE" << endl;
        }
        else {
            resetNodes();
            computePaths(graph.nodes[start]);
            vector<Node*> result = getShortestPathTo(graph.nodes[end]);
            
            printPath(result, start, end);
        }
    }
};

// Get the min of three numbers
int getMin(int a, int b, int c) {
    return min(min(a, b), c);
}

// Binary search for index no recursion
int getCeilingIndex(int array[], int tail[], int left, int right, int key) {
    int mid = -1;
    
    while(right - left > 1) {
        // This is the same as saying (right + left)/2 but safer for a large left
        mid = left + (right - left)/2;
        if(array[tail[mid]] >= key) {
            right = mid;
        }
        else {
            left = mid;
        }
    }
    
    return right;
}

vector<int> getMagiSequenceVector(int array[], int size) {
    
    vector<int> outSequence;
    
    int *tailIndices = new int[size];
    int *prevIndices = new int[size];
    int length;
    
    // Set both arrays to the default values at the start
    memset(tailIndices, 0, sizeof(tailIndices[0])*size);
    memset(prevIndices, 0xFF, sizeof(prevIndices[0])*size);
    
    tailIndices[0] = 0;
    prevIndices[0] = -1;
    length = 1;
    
    for(int i = 1; i < size; i++) {
        if(array[i] < array[tailIndices[0]]) {
            // We have encountered the new smallest value
            tailIndices[0] = i;
        }
        else if(array[i] > array[tailIndices[length-1]]) {
            // Our array wants to extend the longest subsequence
            prevIndices[i] = tailIndices[length-1];
            tailIndices[length++] = i;
        }
        else {
            // Our array element is a possibly in a future subsequence
            int pos = getCeilingIndex(array, tailIndices, -1, length-1, array[i]);
            
            prevIndices[i] = tailIndices[pos-1];
            tailIndices[pos] = i;
        }
    }
    
    // We go through our original array based on the indices from tailIndices
    for(int i = tailIndices[length-1]; i >= 0; i = prevIndices[i]) {
        outSequence.push_back(array[i]);
    }
    
    delete[] tailIndices;
    delete[] prevIndices;
    
    // We will need to flip it back around because it pushes them by backtracking
    std::reverse(outSequence.begin(), outSequence.end());
    
    return outSequence;
}

// This is the edit distance dynamically O(l1 x l2)
int getEditDistance(string word1, int l1, string word2, int l2) {
    
    // Table to store all the results
    int resultsArray[l1 + 1][l2 + 1];
    
    // Fill resultsArray[][] bottom up
    for (int i = 0; i <= l1; i++) {
        for (int j = 0; j <= l2; j++) {
            
            // word1 = "" we have to add all characters from word2 to word1
            if (i == 0) {
                // j = number of operations
                resultsArray[i][j] = j;
            }
            
            // word2 = "" we have to remove all characters from word1
            else if (j == 0) {
                // i = number of operations
                resultsArray[i][j] = i;
            }
            
            // If the last character of each word is the same then move down the words until different
            else if (word1[i-1] == word2[j-1]) {
                resultsArray[i][j] = resultsArray[i-1][j-1];
            }
            
            // If the last character of each word is not the same compute the min frome each operation below
            else {
                // Insert, Remove, Replace
                resultsArray[i][j] = 1 + getMin(resultsArray[i][j-1], resultsArray[i-1][j], resultsArray[i-1][j-1]);
            }
        }
    }
    
    return resultsArray[l1][l2];
}

int main() {
    int numOfRealms = 0;
    cin >> numOfRealms;
    
    Graph graph(numOfRealms);
    
    for (int i = 0; i < numOfRealms; i++) {
        string charmOfRealm = "";
        cin >> charmOfRealm;
        int numOfMagi = 0;
        cin >> numOfMagi;
        int magiArray[numOfMagi];
        
        for (int j = 0; j < numOfMagi; j++) {
            int powerOfMagi = 0;
            cin >> powerOfMagi;
            magiArray[j] = powerOfMagi;
        }
        
        // This will get our longest increasing subsequence
        vector<int> lis = getMagiSequenceVector(magiArray, numOfMagi);
        graph.addNodeToArray(charmOfRealm, i, lis);
        
    }
    
    // call method inside graph class to connect all nodes and set edges
    graph.setupGraph();
    
    string kaelCharm = "";
    cin >> kaelCharm;
    string destCharm = "";
    cin >> destCharm;
    
    int startIndex = graph.findIDWithCharm(kaelCharm);
    int endIndex = graph.findIDWithCharm(destCharm);
    
    // We need to print the proper values for the way there and the way back
    Dijkstra dj(graph);
    dj.doDijkstra(startIndex, endIndex);
    dj.doDijkstra(endIndex, startIndex);
    
    return 0;
}