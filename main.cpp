#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <string.h>
#include <algorithm>

using namespace std;

// Initialize Classes
class Node;
class Edge;
class LinkedList;
class Graph;
class Dijkstra;

// Initialize Methods
vector<int> longestIncreasingSubVector(int array[], int arraySize);
int editDistanceVal(string word1, int l1, string word2, int l2);

// Edge class that connects two nodes
class Edge {
    Node *source;
    Node *destination;
    int weight;
    
public:
    Edge(Node *startNode, Node *endNode, int w) {
        source = startNode;
        destination = endNode;
        weight = w;
    }
    
    Node* getSource() {
        return source;
    }
    
    Node* getDestination() {
        return destination;
    }
    
    int getWeight() {
        return weight;
    }

};

// Node class which represents a realm on the graph
class Node {
    
public:
    string charm;              // Charm for the realm
    int id;                    // Self-Explanatory
    Node *next;                // Stores pointer to next node in linked list for graph
    Node *previous = NULL;     // Stores predecessor once Dijkstra's is performed
    vector<int> lis;           // Longest increasing subsequence for gems needed
    double minDistance = numeric_limits<double>::infinity(); // minDistance for Dijkstra's
    vector<Edge*> edges;       // Vector containing every edge connected to this node
    //bool visited = false;
    
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
    
    // Adds an edge object between this node and an inputted one
    void addEdge(Node *node, int weight) {
        Edge *myEdge = new Edge(this, node, weight);
        edges.push_back(myEdge);
    }
    
    // This prints all the edges and their end destinations (for debugging)
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
    void insertNodeToTail(Node *node) {
        Node* p;
        if (head == NULL)
        {
            p = node;
            head = p;
        }
        else if (tail == NULL)
        {
            p = node;
            tail = p;
            head->next = tail;
        } else {
            p = node;
            tail->next = p;
            tail = p;
        }
    }
    
//    // Prints the charm of each node in the linkedlist (for debugging)
//    void printListOfValues() {
//        Node *p = head;
//        while(p != NULL) {
//            cout << p->charm << " ";
//            p = p->next;
//        }
//    }
    
};

// This is a basic graph class
class Graph {
    int numOfNodes;
    // Array of linked list head pointers containing all adjacent Nodes to a node set at ID
    // i.e. adjNodeList[0] will return a pointer to the head of a linked list of all the
    //      adjacent nodes to node with ID 0
    LinkedList *adjNodeList;
    
public:
    
    // Dynamic array of node pointers
    // It contains all nodes in this particular graph
    // Used an array instead of a matrix or linked list
    // since the graph is not dynamic (extra nodes won't be added once graph is made)
    // Uses ID as index to access nodes; this improves time complexity
    Node* *nodes;

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
        adjNodeList[v1->id].insertNodeToTail(v2);
    }
    
    void addEdgeBidirectional(Node *v1, Node *v2, int weight) {
        v1->addEdge(v2, weight);
        adjNodeList[v1->id].insertNodeToTail(v2);
        
        v2->addEdge(v1, weight);
        adjNodeList[v2->id].insertNodeToTail(v1);
        
    }
    
    // This adds edges between all nodes
    // Time complexity: O(n^2) -> n is number of Nodes in Graph
    void setupGraph() {
        for(int i = 0; i < numOfNodes; i++) {
            for(int j = i + 1; j < numOfNodes; j++) {
                addEdgeBidirectional(nodes[i], nodes[j], editDistanceVal(nodes[i]->charm, (int)nodes[i]->charm.length(), nodes[j]->charm, (int)nodes[j]->charm.length()));
                
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
    
//    // Prints all nodes that can be reached from all nodes (for debugging)
//    void printAllNodes() {
//        for (int i = 0; i < numOfNodes; i++) {
//            cout << "From node " << i << ":\n";
//            adjNodeList[i].printListOfValues();
//            cout << "\n";
//        }
//    }
    
//    // Prints all nodes that can be reached from one node (for debugging)
//    void printNodeFrom(int id) {
//        cout << "From node " << id << ":\n";
//        adjNodeList[id].printListOfValues();
//    }
//    
//    // This prints the nodes a node can see and the weight of each edge between them (for debugging)
//    void printEdges() {
//        for (int i = 0; i < numOfNodes; i++) {
//            cout << "From node " << i << ":\n";
//            nodes[i]->printEdges();
//            cout << "\n";
//        }
//    }
    
//    // This prints the lis of each node (for debugging)
//    void printNodeLis(int id) {
//        vector<int> lis = nodes[id]->lis;
//        
//        for (int i = 0; i < lis.size(); i++) {
//            cout << lis[i] << " ";
//        }
//        cout << "\n";
//    }
};

// Class for Dijkstra's algorithm
class Dijkstra {
    
public:
    
    // We need a reference to the same graph
    Graph graph = NULL;
    
    Dijkstra(Graph g) {
        graph = g;
    }
    
    // This does the Dijkstra algorithm to find the minimum cost path
    // Time complexity: O(|E||V|log(|V|))
    void dijkstraToAllVertices(Node *source) {
        
        // We want a priority queue to have the min at the top
        source->minDistance = 0;
        priority_queue<Node*> nodeQueue;
        
        // We want to push the starting node
        nodeQueue.push(source);
        
        // We go until we reach every node in the single connected graph
        while (!nodeQueue.empty()) {
            Node *u = nodeQueue.top();
            nodeQueue.pop();
            //u->visited = true;
            // Visit each edge that u can reach
            for (Edge *e : u->edges) {
                Node *v = e->getDestination();
                double weight = e->getWeight();
                double distThroughU = u->minDistance + weight;
                
                // If we get a lower distance and also have enough gems to do the incantation
                // This is the greedy step in the algorithm
                // Only changes the min distance if in this iteration is less than the previous
                // and we have enough gems
                if (distThroughU < v->minDistance && u->lis.size() >= e->getWeight()) {
                    
                    priority_queue<Node*> tempQueue;
                    bool done = false;
                    
                    // v is removed from the queue
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
                    v->minDistance = distThroughU;
                    v->previous = u; //predecesor
                    nodeQueue.push(v);
                }
            }
        }
    }
    
    // Chooses a path from the ones returned by the method above
    // Starts at the destination and follows path to the source
    // Returns vector of node pointers in order from source to destination
    vector<Node*> getShortestPathTo(Node *destination) {
        vector<Node*> path;
        for (Node *vertex = destination; vertex != NULL; vertex = vertex->previous) {
            path.push_back(vertex);
        }
        
        // Flips vector around to get the path from the beginning and not from the end
        reverse(path.begin(), path.end());
        
        return path;
    }
    
    //Prints the path from method above
    void printPath(vector<Node*> vec, int start, int end) {
        
        if (graph.nodes[start] != graph.nodes[end] && vec.size() == 1 && vec[vec.size()-1] == graph.nodes[end] && vec[vec.size()-1]->previous == NULL) {
            cout << "IMPOSSIBLE" << endl;
        }
        else {
            int edgeSum = 0;
            int gemSum = 0;
            
            for (int i = 0; i < vec.size() - 1; i++) {
                
                // We get the edit distance between the two nodes
                int editD = editDistanceVal(vec[i]->charm, (int)vec[i]->charm.size(), vec[i + 1]->charm, (int)vec[i + 1]->charm.size());
                edgeSum = edgeSum + editD;
                
                // We sum only the editD number from the longest increasing sequence
                for (int j = 0; j < editD; j++) {
                    gemSum = gemSum + vec[i]->lis[j];
                }
            }
            
            cout << edgeSum << " " << gemSum << endl;
        }
    }
    
    // Resets the nodes to having infinity minDistance and no previous node pointer so
    // we can do Dijkstra's once again
    void resetNodes() {
        for (int i = 0; i < graph.pubNumOfNodes; i++) {
            graph.nodes[i]->previous = NULL;
            graph.nodes[i]->minDistance = numeric_limits<double>::infinity();
            //graph.nodes[i]->visited = false;
        }
    }
    
    // Method to put all of the above methods together
    void performDijstras(int start, int end) {
        if (start == -1 || end == -1) {
            cout << "IMPOSSIBLE" << endl;
        }
        else {
            resetNodes();
            dijkstraToAllVertices(graph.nodes[start]);
            vector<Node*> result = getShortestPathTo(graph.nodes[end]);
            
            printPath(result, start, end);
        }
    }
};

// Get the min of three numbers
int getMin(int a, int b, int c) {
    return min(min(a, b), c);
}

// Binary search for index
int getCeilingIndex(int array[], int tail[], int left, int right, int key) {
    int mid = -1;
    
    while(right - left > 1) {
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

// Returns vector with longest increasing subsequence of inputted int array
// Time Complexity: O(n^2) -> n is size of input array
vector<int> longestIncreasingSubVector(int array[], int size) {
    
    vector<int> outputSequence;
    
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
        outputSequence.push_back(array[i]);
    }
    
    delete[] tailIndices;
    delete[] prevIndices;
    
    // Flips vector back around because it pushes them into the vector by backtracking
    reverse(outputSequence.begin(), outputSequence.end());
    
    return outputSequence;
}

// Performs edit distance dynamically
// Time complexity: O(l1 x l2)
int editDistanceVal(string word1, int l1, string word2, int l2) {
    
    // Table to store all the results
    int resultsArray[l1 + 1][l2 + 1];
    
    // Fill resultsArray[][] bottom up
    for (int i = 0; i <= l1; i++) {
        for (int j = 0; j <= l2; j++) {
            
            // if word1 = "" (empty) we have to add all characters from word2 to word1
            if (i == 0) {
                // j = number of operations
                resultsArray[i][j] = j;
            }
            
            // if word2 = "" (empty) we have to remove all characters from word1
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
                //                              Insert,               Remove,               Replace
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
        
        // Longest increasing subsequence of magi from inputted realm
        vector<int> lis = longestIncreasingSubVector(magiArray, numOfMagi);
        graph.addNodeToArray(charmOfRealm, i, lis);
        
    }
    
    // Connects all nodes and makes edges
    graph.setupGraph();
    //graph.printEdges();
    //graph.printAllNodes();
    
    string kaelCharm = "";
    cin >> kaelCharm;
    string destCharm = "";
    cin >> destCharm;
    
    int startIndex = graph.findIDWithCharm(kaelCharm);
    int endIndex = graph.findIDWithCharm(destCharm);
    
    // We need to print the proper values for the way there and the way back
    Dijkstra dij(graph);
    dij.performDijstras(startIndex, endIndex);
    dij.performDijstras(endIndex, startIndex);
    
    return 0;
}
