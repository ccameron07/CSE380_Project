#ifndef ELEMENTS_H
#define ELEMENTS_H
#include <vector>
#include "./eigen3/Eigen/Dense"

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

class Node1d : public Node {
    public:
        double coords;
        //Constructor Method
        Node1d(int n_init, double x, int df_init = 0, bool boundary_init = false, double BC_init = 0.0);
 };

class Node2d : Node {
    public:
        Vector2d coords ;
        //Constructor Method
        Node2d(int n_init, double x, double y, int df_init = 0, bool boundary_init = false, double BC_init = 0.0) ;
};

class Line {
    public:
        vector<Node*> nodes;
        
        Line(Node& n1, Node& n2) ;

        void addNodes(int n, vector<Node>& nodes_g, vector<Node*>& nodes_e) ;
};

#endif
