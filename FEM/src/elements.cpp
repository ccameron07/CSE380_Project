#include <vector>
#include "./eigen3/Eigen/Dense"
#include "elements.hpp"

/*
using namespace Eigen ;
using namespace std ;
class Node {
    private:
        Node(){ } 
    public:
        int n, df ;
        bool boundary ;
        double BC ;
        Node ( int n_init, int df_init = 0, bool boundary_init = false, double BC_init = 0.0 );
};
*/

Node::Node(int n_init, int df_init, bool boundary_init, double BC_init) {
    n = n_init ;
    df = df_init ;
    boundary = boundary_init;
    BC = BC_init ;
}

Node1d::Node1d(int n_init, double x, int df_init, bool boundary_init, double BC_init) 
:Node(n_init, df_init, boundary_init, BC_init)
{
    coords = x ;
}

Node2d::Node2d(int n_init, double x, double y, int df_init, bool boundary_init, double BC_init) 
:Node(n_init, df_init, boundary_init, BC_init)
{
    coords(0) = x ;
    coords(1) = y ;
}

Line::Line(Node* n1, Node* n2) {
            nodes.push_back(n1);
            nodes.push_back(n2);
}

void Line::addNodes(int n, vector<Node>& nodes_g, vector<Node*>& nodes_e) {
     
     if (nodes.size() == (2+n)) {
     
         for(int i = nodes.size(); i>2; i--){
             nodes_e.push_back(nodes[i-1]) ;
         }

     } else {
         int n_nodes = nodes_g.size() ;
         
         for(int i = 0; i < n; i++){
             nodes_g.push_back(Node(n_nodes+i, 0, (*nodes[0]).boundary, (*nodes[0]).BC)) ;
             nodes_e.push_back(&nodes_g[n_nodes]) ;
             nodes.push_back(&nodes_g[n_nodes]) ;
         }
     }
 }
