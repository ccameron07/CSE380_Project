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
