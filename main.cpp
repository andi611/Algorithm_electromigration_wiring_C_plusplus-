/****************************************************************************
  FileName     [ main.cpp ]
  Synopsis     [ Implementation of a Optimal Wiring Topology EM Router ]
  Author       [ Ting-Wei (Andy) Liu ]
  Copyright    [ Copyleft(c), NTUEE, NTU, Taiwan ]
****************************************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdio>
#include <ctype.h>
#include <cassert>
#include <cstring>
#include <climits>
#include <cmath>
#include <map>
#include <float.h>
#include <utility>
#include <algorithm>
using namespace std;

/**********************************/
/*          Declarations          */
/**********************************/
class Edge;
class Node;
typedef pair<int, int>  coordinate;
static string str_buf;
#define GREEDY
//define DEBUG


/****************************/
/*        class Edge        */
/****************************/
class Edge
{
public:
    Edge() : _flow(INT_MAX), _original_edge(NULL) {}
    

    // access functions
    Node* get_node_u() const { return _node_u; }
    Node* get_node_v() const { return _node_v; }
    int get_wirelength() const { return _wirelength; }
    int get_flow() const { return _flow; }
    int get_capacity() const { return _capacity; }
    Edge* get_original_edge() const { return _original_edge; }

    // set functions
    void set_node(Node* u, Node* v) { _node_u = u; _node_v = v; }
    void set_node_u(Node* u) { _node_u = u; }
    void set_node_v(Node* v) { _node_v = v; }
    void set_wirelength(int w) { _wirelength = w; }
    void set_flow(int f) { _flow = f; }
    void set_capacity(int c) { _capacity = c; }
    void set_original_edge(Edge* e) { _original_edge = e; }

    // set all function
    void set_edge(Node* u, Node* v, int w, int f, int c) { 
            _node_u = u; 
            _node_v = v; 
            _wirelength = w;
            _flow = f;
            _capacity = c; }

private:
    Node*           _node_u; // directed edge from u to v
    Node*           _node_v;
    int             _wirelength;
    int             _flow;
    int             _capacity;
    Edge*           _original_edge; // match residual edges to original flow edges
};


/****************************/
/*        class Node        */
/****************************/
class Node
{
public:
    Node(coordinate c, int f) : _coordinate(c), _max_flow(f), _cur_flow(0) {}
    
    // basic operations
    void set_index(int i) { _idx = i; }
    int get_index() const { return _idx; }
    void set_max_flow(int f) { _max_flow = f; }
    void set_cur_flow(int f) { _cur_flow = f; }
    int get_max_flow() const { return _max_flow; }
    int get_cur_flow() const { return _cur_flow; }
    coordinate get_coordinate() const { return _coordinate; }
    int check_flow_left() const { return _max_flow - _cur_flow; }
    
    // fanout operations
    bool fanout_empty() const { return _fanout.empty(); }
    size_t get_fanout_size() const { return _fanout.size(); }

    // access fanout
    Edge* get_fanout(size_t idx) const { return _fanout.at(idx); }
    void add_fanout(Edge* e) { _fanout.push_back(e); }
    void clear_fanout() { _fanout.clear(); }

private: 
    // node information          
    int             _idx;
    coordinate      _coordinate;
    int             _max_flow;
    int             _cur_flow;
    // node connections
    vector<Edge*>   _fanout;
};


/*********************************/
/*        class EM Router        */
/*********************************/
class EM_Router
{
public:
    EM_Router() {
        _dummy_s = new Node(coordinate(0,0), 0);
        _dummy_t = new Node(coordinate(0,0), 0);
        _dummy_s_res = new Node(coordinate(0,0), 0);
        _dummy_t_res = new Node(coordinate(0,0), 0);
        _T = new Node(coordinate(0,0), 0);
        _T->set_index(0);
    }

    // mian functions for flow net
    int build_flow_network();
    void connect_flow_network();
    bool push_flow();
    bool push_flow_greedy();

    // main functions for residual net
    void build_residual_network();
    void connect_residual_network();
    void disconnect_residual_network();

    // neg cycle operations
    vector<Edge*> find_neg_cycle();
    void augment_graph_addT(Node* T);
    void push_flow_neg_cycle(vector<Edge*>& neg_cycle);
    void preprocess_node();

    // utility functions
    bool read_input_table (const string&);
    void display_input (bool show_input = false) const;
    void display_graph (bool show_graph = false) const;
    int calculate_wire_area() const;
    bool check_valid_cycle(vector<Edge*> neg_cycle) const;
    bool check_valid_flow(int total_flow) const;
    void save_output(const string& output_file) const;

private:
    // flow network
    int                     _n_source_and_sink;
    vector< vector<int> >   _input_table;
    Node*                   _dummy_s;
    Node*                   _dummy_t;
    // residual network
    Node*                   _dummy_s_res;
    Node*                   _dummy_t_res;
    Node*                   _T;
    vector<Node*>           _node_list;
    int                     _total_nodes;

    // helper functions
    int calculate_wire_len(int x1, int x2, int y1, int y2) { return (abs(x1 - x2) + abs(y1 - y2)); }
    int min_on_int(int int1, int int2) { if (int1 <= int2) return int1; else return int2; }
    int is_in_vector_idx(int v, vector<int>& vec);
 };


/*******************************/
/*        class Compare        */
/*******************************/
class Compare_Edge { public: bool operator() (Edge* e1, Edge* e2) { return (e1->get_wirelength() < e2->get_wirelength()); } };
class Compare_Vec { 
    public: bool operator() (vector<int> v1, vector<int> v2) {
        if (v1.at(0) < v2.at(0)) return true;
        if (v1.at(0) > v2.at(0)) return false;
        if (v1.at(1) < v2.at(1)) return true;
        if (v1.at(1) > v2.at(1)) return false;
        if (v1.at(2) < v2.at(2)) return true;
        if (v1.at(2) > v2.at(2)) return false;
        if (v1.at(3) < v2.at(3)) return true;
        if (v1.at(3) > v2.at(3)) return false;
        if (v1.at(4) < v2.at(4)) return true;
        if (v1.at(4) > v2.at(4)) return false;
        return false;
        }
};


/************************************************/
/*        class EM Router main functions        */
/************************************************/
int EM_Router::build_flow_network ()
{      
    int n_source = 0;
    int n_sink = 0;
    // build nodes
    for (int i = 0; i < _input_table.size(); ++i) {
        if (_input_table.at(i).at(2) > 0) {
            // initiate node
            coordinate xy(_input_table.at(i).at(0), _input_table.at(i).at(1));
            Node* node_s = new Node(xy, abs(_input_table.at(i).at(2)));

            // set edges
            Edge* edge = new Edge();
            edge->set_node(_dummy_s, node_s);

            // connect to dummy node
            _dummy_s->add_fanout(edge);
            n_source += 1;
        }
        else {
            // initiate node
            coordinate xy(_input_table.at(i).at(0), _input_table.at(i).at(1));
            Node* node_t = new Node(xy, abs(_input_table.at(i).at(2)));

            // set edges
            Edge* edge = new Edge();
            edge->set_node(_dummy_t, node_t);

            // connect to dummy node
            _dummy_t->add_fanout(edge);
            n_sink += 1;
        }
    }
    // connect s and t nodes
    connect_flow_network();
    _total_nodes = n_source + n_sink + 1;

    // calculate total possible for dummy nodes
    int total_flow_s = 0;
    int total_flow_t = 0;
    for (int i = 0; i < _dummy_s->get_fanout_size(); ++i) { total_flow_s += _dummy_s->get_fanout(i)->get_node_v()->get_max_flow(); }
    for (int i = 0; i < _dummy_t->get_fanout_size(); ++i) { total_flow_t += _dummy_t->get_fanout(i)->get_node_v()->get_max_flow(); }
    _dummy_s->set_max_flow(total_flow_s);
    _dummy_t->set_max_flow(total_flow_t);

    #ifdef DEBUG
    assert (total_flow_s == total_flow_t);
    #endif
    
    cout << "Number of source: " << n_source << endl;
    cout << "Number of sink: " << n_sink << endl;
    cout << "Complete building flow network with total flow: " << total_flow_s << endl;
    return total_flow_s;
}


void EM_Router::connect_flow_network ()
{   
    for (int i = 0; i < _dummy_s->get_fanout_size(); ++i) {
        Node* node_s = _dummy_s->get_fanout(i)->get_node_v();
        
        for (int j = 0; j < _dummy_t->get_fanout_size(); ++j) {
            // initiate
            Node* node_t = _dummy_t->get_fanout(j)->get_node_v();
            Edge* edge = new Edge();

            // calculate wirelength
            coordinate xy_s =  node_s->get_coordinate();
            coordinate xy_t =  node_t->get_coordinate();
            int len = calculate_wire_len(xy_s.first, xy_t.first, xy_s.second, xy_t.second);

            // calculate capacity
            int flow_s = node_s->get_max_flow();
            int flow_t = node_t->get_max_flow();
            int cap = min_on_int(flow_s, flow_t);

            // set edges
            edge->set_edge(node_s, node_t, len, INT_MAX, cap);

            // connect
            node_s->add_fanout(edge);
        }
    }
}


bool EM_Router::push_flow_greedy () 
{   
    bool pushable = false;
    cout << "Pushing initial flow with greedy Method!" << endl;
    // colllect all edges
    vector<Edge*> edge_list;
    for (int i = 0; i < _dummy_s->get_fanout_size(); ++i) {
        Node* node_s = _dummy_s->get_fanout(i)->get_node_v();
        for (int j = 0; j < node_s->get_fanout_size(); ++j)       
            edge_list.push_back(node_s->get_fanout(j));
    }

    std::sort(edge_list.begin(), edge_list.end(), Compare_Edge());

    // greedy push: starting with smallest wire length
    Node* node_s;
    Node* node_t;
    int flow_s;
    int flow_t;
    for (int i = 0; i < edge_list.size(); ++i) {
        node_s = edge_list.at(i)->get_node_u();
        node_t = edge_list.at(i)->get_node_v();
        flow_s = node_s->check_flow_left();
        if (flow_s > 0) {   
            flow_t = node_t->check_flow_left();
            if (flow_t > 0) {
                pushable = true;
                int push_flow = min_on_int(flow_s, flow_t);

                // push flow
                node_s->set_cur_flow(node_s->get_cur_flow()+push_flow);
                node_t->set_cur_flow(node_t->get_cur_flow()+push_flow);
                edge_list.at(i)->set_flow(push_flow);

                #ifdef DEBUG
                assert(node_s->get_cur_flow() <= node_s->get_max_flow());
                assert(node_t->get_cur_flow() <= node_t->get_max_flow());
                #endif
                }
            }
    }

    return pushable;
}



bool EM_Router::push_flow () 
{   
    #ifdef GREEDY
    if (_total_nodes > 100) return push_flow_greedy();
    #endif
    cout << "Pushing initial flow with normal Method." << endl;

    bool pushable = false;
    for (int i = 0; i < _dummy_s->get_fanout_size(); ++i) {
        // node s
        Node* node_s = _dummy_s->get_fanout(i)->get_node_v();

        for (int j = 0; j < node_s->get_fanout_size(); ++j) {
            
            Edge* edge = node_s->get_fanout(j);
            int flow_s = node_s->check_flow_left();
            if (flow_s > 0) {   
                // nove t
                Node* node_t = edge->get_node_v();

                // calculate flow
                int flow_t = node_t->check_flow_left();
                if (flow_t > 0) {
                    pushable = true;
                    int push_flow = min_on_int(flow_s, flow_t);

                    // push flow
                    node_s->set_cur_flow(node_s->get_cur_flow()+push_flow);
                    node_t->set_cur_flow(node_t->get_cur_flow()+push_flow);
                    edge->set_flow(push_flow);

                    #ifdef DEBUG
                    assert(node_s->get_cur_flow() <= node_s->get_max_flow());
                    assert(node_t->get_cur_flow() <= node_t->get_max_flow());
                    #endif
                }
            }      
        }  
    }
    return pushable;
}


void EM_Router::build_residual_network ()
{   
    cout << "Building residual network, starting negative cycle removal..." << endl;
    int idx = 1;
    for (int i = 0; i < _dummy_s->get_fanout_size(); ++i) {
        // get origninal node and assign idx
        Node* node_o = _dummy_s->get_fanout(i)->get_node_v();
        node_o->set_index(idx);
        idx += 1;
        
        // copy node
        Node* node_s = new Node(*node_o); 
        node_s->clear_fanout();

        // set edge
        Edge* edge = new Edge();
        edge->set_node(_dummy_s_res, node_s);

        // connect to dummy node
        _dummy_s_res->add_fanout(edge);
    }     
    for (int i = 0; i < _dummy_t->get_fanout_size(); ++i) {
        // get origninal node and assign idx
        Node* node_o = _dummy_t->get_fanout(i)->get_node_v();
        node_o->set_index(idx);
        idx += 1;

        // copy node
        Node* node_t = new Node(*node_o); 
        node_t->clear_fanout();

        // set edge 
        Edge* edge = new Edge();
        edge->set_node(_dummy_t_res, node_t);

        // connect to dummy node
        _dummy_t_res->add_fanout(edge);
    }
}


void EM_Router::connect_residual_network()
{   
    map<Edge*, Edge*> edge2edge; // map residual edge to original edge

    for (int i = 0; i < _dummy_s->get_fanout_size(); ++i){
        Node* node_s = _dummy_s->get_fanout(i)->get_node_v();
        Node* node_s_res = _dummy_s_res->get_fanout(i)->get_node_v();

        for (int j = 0; j < node_s->get_fanout_size(); ++j) {
            Edge* edge_o = node_s->get_fanout(j);
            Node* node_t_res = _dummy_t_res->get_fanout(j)->get_node_v();
            
            #ifdef DEBUG
            assert (edge_o->get_wirelength() != INT_MAX);
            #endif

            // forward edge
            if ((edge_o->get_flow() < edge_o->get_capacity()) || (edge_o->get_flow() == INT_MAX)) {
                Edge* edge_f = new Edge();
                edge_f->set_edge(node_s_res, node_t_res, edge_o->get_wirelength(), INT_MAX, edge_o->get_capacity()-edge_o->get_flow());
                edge_f->set_original_edge(edge_o);
                node_s_res->add_fanout(edge_f);
            }
            // backward edge 
            if ((edge_o->get_flow() > 0) && (edge_o->get_flow() != INT_MAX)) {
                Edge* edge_b = new Edge();
                edge_b->set_edge(node_t_res, node_s_res, -edge_o->get_wirelength(), INT_MAX, -edge_o->get_flow());
                edge_b->set_original_edge(edge_o);
                node_t_res->add_fanout(edge_b);
            }
        }
    }
}


void EM_Router::disconnect_residual_network()
{
    for (int i = 0; i < _dummy_s_res->get_fanout_size(); ++i){
        Node* node_s = _dummy_s_res->get_fanout(i)->get_node_v();
        for (int j = 0; j < node_s->get_fanout_size(); ++j){
            delete node_s->get_fanout(j); // delete edge
        }
        node_s->clear_fanout(); // delete fanout
    }
    for (int i = 0; i < _dummy_t_res->get_fanout_size(); ++i){
        Node* node_t = _dummy_t_res->get_fanout(i)->get_node_v();
        for (int j = 0; j < node_t->get_fanout_size(); ++j){
            delete node_t->get_fanout(j); // delete edge
        }
        node_t->clear_fanout(); // delete fanout
    }
}


/******************************************************/
/*        class EM Router neg cycle operations        */
/******************************************************/
vector<Edge*> EM_Router::find_neg_cycle()
{
    // implementation of the Bellman-Ford Dynamic Programming algorithm

    // Augment Graph, add a new node T, connect all nodes to T
    augment_graph_addT(_T);

    // construct memories and trace backs
    int* memory_opt = new int[_total_nodes];
    int* memory_old = new int[_total_nodes];
    int* trace_node = new int[_total_nodes];
    Edge** trace_edge = new Edge*[_total_nodes];

    for(int i = 0; i < _total_nodes; ++i) {
        memory_opt[i] = INT_MAX-77777;
        memory_old[i] = INT_MAX-77777;
        trace_node[i] = 0;
        trace_edge[i] = NULL;
    }
    memory_old[0] = 0;

    int* holder;
    int edge_cost;
    int idx;

    for (int i = 1; i < _total_nodes; ++i){
        for (int j = 0; j < _total_nodes; ++j){
            
            // find min cost over all edges of node j
            memory_opt[j] = memory_old[j];
            for (int k = 0; k < _node_list.at(j)->get_fanout_size(); ++k) {
                idx = _node_list.at(j)->get_fanout(k)->get_node_v()->get_index();
                edge_cost = _node_list.at(j)->get_fanout(k)->get_wirelength() + memory_old[idx];
                if (edge_cost < memory_opt[j]) {
                    // record optimal
                    memory_opt[j] = edge_cost;
                    // record trace back
                    trace_node[j] = idx;
                    trace_edge[j] = _node_list.at(j)->get_fanout(k);
                }
            }
        }
        // swap
        holder = memory_old;
        memory_old = memory_opt;
        memory_opt = holder;
    }

    // check if odd cycle exist
    int neg_cycle_idx = INT_MAX;
    for (int i = 0; i < _total_nodes; ++i)
        if (memory_opt[i] != memory_old[i]) {
            neg_cycle_idx = i;
            break;
        }


    vector<Edge*> neg_cycle_edge;
    vector<int> neg_cycle_node;

    // if no neg cycle
    if (neg_cycle_idx == INT_MAX){
        return neg_cycle_edge;
    }

    // if neg cycle
    for (int i = _total_nodes-1; i >= 0; --i){
        neg_cycle_node.push_back(neg_cycle_idx);
        neg_cycle_edge.push_back(trace_edge[neg_cycle_idx]);
        neg_cycle_idx = trace_node[neg_cycle_idx];
        idx = is_in_vector_idx(neg_cycle_idx, neg_cycle_node);
        if (idx != INT_MAX) break;
    }
    // truncate the front
    neg_cycle_edge.erase(neg_cycle_edge.begin(),neg_cycle_edge.begin()+idx);
    
    delete[] memory_opt;
    delete[] memory_old;
    delete[] trace_node;
    delete[] trace_edge;
    return neg_cycle_edge;
}


void EM_Router::push_flow_neg_cycle(vector<Edge*>& neg_cycle)
{   
    int push_flow = INT_MAX;
    for (int i = 0; i < neg_cycle.size(); ++i){
        if (abs(neg_cycle.at(i)->get_capacity()) < push_flow)
            push_flow = abs(neg_cycle.at(i)->get_capacity());
    }
    Edge* edge_o = NULL;
    for (int i = 0; i < neg_cycle.size(); ++i){
        edge_o = neg_cycle.at(i)->get_original_edge();
        
        if (edge_o->get_flow() == INT_MAX) edge_o->set_flow(0); // initial flow if unset
        if (neg_cycle.at(i)->get_wirelength() < 0) // a backward edge
            edge_o->set_flow(edge_o->get_flow()-push_flow);
        else
            edge_o->set_flow(edge_o->get_flow()+push_flow);
        
        #ifdef DEBUG
        assert(edge_o->get_flow() <= edge_o->get_capacity());
        assert(edge_o->get_flow() >= 0);
        #endif
    }
}



/**************************************************/
/*        class EM Router helper functions        */
/**************************************************/
void EM_Router::augment_graph_addT(Node* T)
{
    for (int i = 0; i < _dummy_s_res->get_fanout_size(); ++i){
        Node* node_s = _dummy_s_res->get_fanout(i)->get_node_v();
        Edge* edge_T = new Edge();
        edge_T->set_edge(node_s, T, 0, INT_MAX, 0); // 0 cost edge
        node_s->add_fanout(edge_T);
    }
    for (int i = 0; i < _dummy_t_res->get_fanout_size(); ++i){
        Node* node_t = _dummy_t_res->get_fanout(i)->get_node_v();
        Edge* edge_T = new Edge();
        edge_T->set_edge(node_t, T, 0, INT_MAX, 0); // 0 cost edge
        node_t->add_fanout(edge_T);
    }
}


void EM_Router::preprocess_node()
{
    // construct node list
    _node_list.push_back(_T);
    for (int i = 0; i < _dummy_s_res->get_fanout_size(); ++i)
        _node_list.push_back(_dummy_s_res->get_fanout(i)->get_node_v());
    for (int i = 0; i < _dummy_t_res->get_fanout_size(); ++i)
        _node_list.push_back(_dummy_t_res->get_fanout(i)->get_node_v());
    #ifdef DEBUG
    assert(_total_nodes == _node_list.size());
    #endif
}


int EM_Router::is_in_vector_idx(int v, vector<int>& vec)
{
    for (int i = 0; i < vec.size(); ++i)
        if (vec.at(i) == v)
            return i;
    return INT_MAX;
}


/***************************************************/
/*        class EM Router utility functions        */
/***************************************************/
bool EM_Router::read_input_table (const string& input_file)
{
    // opening file
    ifstream ifs(input_file.c_str());
    if (!ifs) {
        cerr << "Cannot open file \"" << input_file << "\"!" << endl;
        return false; }

    // get first line
    getline(ifs, str_buf, '\n');
    _n_source_and_sink = atoi(str_buf.c_str());

    // get middle line
    for (int i=0; i<_n_source_and_sink; ++i) {
        getline(ifs, str_buf, '\n');
        
        int begin = 0;
        int end = 1;
        vector<int> line;
        
        // parse each token
        for (string::iterator it=str_buf.begin(); it!=str_buf.end(); ++it) {
            if ((*it == ' ') || (*it == '\t')) {
                line.push_back(atoi(str_buf.substr(begin, end - begin).c_str()));
                begin = end;
            }
            end += 1;
        }

        // get the last token
        line.push_back(atoi(str_buf.substr(begin, end - begin).c_str()));
        _input_table.push_back(line);
    }
    return true;
}


void EM_Router::display_input(bool show_input) const
{
    if (show_input == true) {
        for (int i = 0; i < _input_table.size(); ++i) {
            for (int j = 0; j < _input_table.at(i).size(); ++j) {
                cout << _input_table.at(i).at(j) << ' ';
            }
            cout << endl;
        }
    }
    #ifdef DEBUG
    assert(_n_source_and_sink == _input_table.size());
    #endif
    cout << "Number of sources + sinks: " << _n_source_and_sink << endl;
}


void EM_Router::display_graph(bool show_graph) const
{
    if (show_graph == true) {
        for (int i = 0; i < _dummy_s->get_fanout_size(); ++i){
            Node* node_s = _dummy_s->get_fanout(i)->get_node_v();

            for (int j = 0; j < node_s->get_fanout_size(); ++j) {
                Edge* edge = node_s->get_fanout(j);
                cout << "S: " << i+1 << " T: " << j+1 << "      Flow: " << edge->get_flow() << " / Capacity: " << edge->get_capacity() << endl;
            }
        }
    }
}


int EM_Router::calculate_wire_area() const
{   
    int area = 0;
    for (int i = 0; i < _dummy_s->get_fanout_size(); ++i) {
        Node* node_s = _dummy_s->get_fanout(i)->get_node_v();
        for (int j = 0; j < node_s->get_fanout_size(); ++j) {
            if (node_s->get_fanout(j)->get_flow() != INT_MAX)
                area += node_s->get_fanout(j)->get_wirelength() * node_s->get_fanout(j)->get_flow();
        }
    }
    cout << "Final Area: " << area << endl;
    return area;
}


bool EM_Router::check_valid_cycle(vector<Edge*> neg_cycle) const
{
    if (neg_cycle.size() == 0) return true;
    int sum = 0;
    for (int i = 0; i < neg_cycle.size(); ++i) {
        sum += neg_cycle.at(i)->get_wirelength();
    }
    if (sum < 0) return true;
    else return false;
}


bool EM_Router::check_valid_flow(int total_flow) const
{   
    // fill in flow table
    int idx_s;
    int idx_t;
    vector< vector<int> > flow = vector< vector<int> >(_dummy_s->get_fanout_size(), vector<int>(_dummy_t->get_fanout_size(), 0));
    for (int i = 0; i < _dummy_s->get_fanout_size(); ++i) {
        Node* node_s = _dummy_s->get_fanout(i)->get_node_v();
        for (int j = 0; j < node_s->get_fanout_size(); ++j) {
            if (node_s->get_fanout(j)->get_flow() != INT_MAX) {
                idx_s = node_s->get_index() - 1;
                idx_t = node_s->get_fanout(j)->get_node_v()->get_index() - _dummy_s->get_fanout_size() - 1;
                flow.at(idx_s).at(idx_t) = node_s->get_fanout(j)->get_flow();
            }
        }
    }

    int idx = 0;
    int sum_s = 0;
    int sum_t = 0;
    // check the sum of flow of source
    for (int i = 0; i < _dummy_s->get_fanout_size(); ++i) {
        for (int j = 0; j < _dummy_t->get_fanout_size(); ++j) {
            sum_s += flow.at(i).at(j);
        }
        idx += 1;
        if (sum_s != _node_list.at(idx)->get_max_flow()) return false;
        sum_s = 0;
    }
    // check the sum of flow of sink
    for (int i = 0; i < _dummy_t->get_fanout_size(); ++i) {
        for (int j = 0; j < _dummy_s->get_fanout_size(); ++j) {
            sum_t += flow.at(j).at(i);
        }
        idx += 1;
        if (sum_t != _node_list.at(idx)->get_max_flow()) return false;
        sum_t = 0;
    }
    return true;
}


void EM_Router::save_output(const string& output_file) const {
    ofstream file;
    file.open (output_file.c_str());
    int wire_area = calculate_wire_area();
    file << wire_area << '\n';

    vector<vector<int> > result;
    for (int i = 0; i < _dummy_s->get_fanout_size(); ++i) {
        Node* node_s = _dummy_s->get_fanout(i)->get_node_v();
        for (int j = 0; j < node_s->get_fanout_size(); ++j) {
            Edge* edge = node_s->get_fanout(j);
            if ((edge->get_flow() != INT_MAX) && (edge->get_flow() > 0)) {
                Node* node_u = edge->get_node_u();
                Node* node_v = edge->get_node_v();
                vector<int> row = vector<int>(5, INT_MAX);

                row.at(0) = node_u->get_coordinate().first;
                row.at(1) = node_u->get_coordinate().second;
                row.at(2) = node_v->get_coordinate().first;
                row.at(3) = node_v->get_coordinate().second;
                row.at(4) = edge->get_flow();
                result.push_back(row);
            }
        }
    }
    std::sort(result.begin(), result.end(), Compare_Vec());

    // output
    for (int i = 0; i < result.size(); ++i) {
        for (int j = 0; j < result.at(0).size(); ++j) {
            if (j == result.at(0).size() - 1) file << result.at(i).at(j);
            else file << result.at(i).at(j) << ' ';
        }
        file << '\n';
    }
    file.close();
    cout << "Result successfully saved to: " << output_file << endl;
}


/**********************************/
/*          Main Function         */
/**********************************/
int main(int argc, char** argv)
{   
    // Settings
    bool show_input = false;
    bool show_graph = false;

    if (argc != 3) { // argument count
        cout << "Illegal number of arguments! Terminating program!" << endl; 
        return 1;
    }
    
    EM_Router* EMR = new EM_Router();

    EMR->read_input_table(argv[1]);
    EMR->display_input(show_input);
    EMR->display_graph(show_graph);

    // flow network
    int total_flow = EMR->build_flow_network();
    // push flow
    EMR->push_flow();
    
    EMR->build_residual_network();
    EMR->preprocess_node();

    int iteration = 0;
    while (true) {
        ++iteration;
        cout << "Iteration: " << iteration << endl;
        
        // connect residual
        EMR->connect_residual_network();

        // find and push neg cycle
        vector<Edge*> neg_cycle = EMR->find_neg_cycle();

        #ifdef DEBUG
        assert(EMR->check_valid_cycle(neg_cycle));
        #endif

        if (neg_cycle.size() > 0)
            EMR->push_flow_neg_cycle(neg_cycle);

        // disconnect residual
        EMR->disconnect_residual_network();
        if (neg_cycle.size() == 0) break;
    }

    #ifdef DEBUG
    assert(EMR->check_valid_flow(total_flow));
    #endif
    EMR->save_output(argv[2]);


    return 0;
}