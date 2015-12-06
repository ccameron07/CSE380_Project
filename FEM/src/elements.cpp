#include <vector>
#include "./eigen3/Eigen/Dense"
#include "elements.hpp"
#include <iostream>
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
    ind = n_init ;
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

Line1d::Line1d(Node1d* n1, Node1d* n2) {
            nodes.push_back(n1);
            nodes.push_back(n2);
}

void Line1d::addNodes(int n, vector<Node1d> &nodes_g, vector<Node1d*> &nodes_e) {
     
     if (nodes.size() == (2+n)) {
     
         for(int i = nodes.size(); i>2; i--){
             nodes_e.push_back(nodes[i-1]) ;
         }

     } else {
         int n_nodes = nodes_g.size() ;
         
         for(int i = 0; i < n; i++){
             nodes_g.push_back(Node1d(n_nodes+i, 0, (*nodes[0]).boundary, (*nodes[0]).BC)) ;
             nodes_e.push_back(&nodes_g[n_nodes]) ;
             nodes.push_back(&nodes_g[n_nodes]) ;
         }
     }
}

Element1d::Element1d(int ind_init, int order_init, int quad_init, vector<Node1d*> Nodes_init, vector<Line1d*> Edges_init){
    ind = ind_init ;
    order = order_init ;
    quad_pts = quad_init ;
    Nodes = Nodes_init ;
    Edges = Edges_init ;
    jacobian_calc() ;
}

void Element1d::jacobian_calc(){
    jac_det = 0.5*(Nodes[1]->coords - Nodes[0]->coords) ;
}

double Element1d::master_2_global( double xi ) {
    double x = Nodes[0]->coords + (1.0 + xi)*jac_det ;
    return x ;
}

void Element1d::Element1d::quadrature(int n) {
    int ind = (quad_pts-1)*quad_pts/2 + n ;
    vector<double> ws {2.0,1.0, 1.0,
                            0.8888888888888888888888889,
                            0.5555555555555555555555556,
                            0.5555555555555555555555556,
                            0.6521451548625461426269361,
                            0.6521451548625461426269361,
                            0.3478548451374538573730639,
                            0.3478548451374538573730639} ;
    
    vector<double> xis {0.0,
                            0.5773502691896257645091488,
                            -0.5773502691896257645091488,
                            0.0,
                            0.7745966692414833770358531,
                            -0.7745966692414833770358531,
                            0.3399810435848562648026658,
                            -0.3399810435848562648026658,
                            0.8611363115940525752239465,
                            -0.8611363115940525752239465} ;
    w = ws[ind];
    xi = xis[ind];
}

double Element1d::shape(int n, double xi) {
    double psi ;
    
    switch(order){
        case 1:
            switch(n){
                case 0:
                    psi = 0.5*(1.0-xi);
                    return psi ;
                    break ;
                case 1:
                    psi = 0.5*(1.0+xi);
                    return psi ;
                    break ;
            }

        case 2:
            switch(n){
                case 0:
                    psi = 0.5*xi*(xi-1.0);
                    return psi ;
                    break ;
                case 1:
                    psi = (1.0-xi*xi);
                    return psi ;
                    break ;
                case 2:
                    psi = 0.5*xi*(xi+1.0);
                    return psi ;
                    break ;           
            }
    }
}

double Element1d::dshape(int n, double xi) {
    double dpsi ;

    switch(this->order){
        case 1:
            switch(n){
                case 0:
                    dpsi = -0.5;
                    return dpsi/jac_det ;
                    break ;
                case 1:
                    dpsi = 0.5;
                    return dpsi/jac_det ;
                    break ;
            }

        case 2:
            switch(n){
                case 0:
                    dpsi = xi-0.5;
                    return dpsi/jac_det ;
                    break ;
                case 1:
                    dpsi = 1.0-2.0*xi;
                    return dpsi/jac_det ;
                    break ;
                case 2:
                    dpsi = xi+0.5;
                    return dpsi/jac_det ;
                    break ;           
            }
    }
}

void Element1d::kfCalc(MatrixXd &A, VectorXd &b) {

    int size = A.rows() ;
    double stiffness = 1.0 ;
    double force = 1.0 ;
    int i_g, j_g ;

    for(int i = 0 ; i < order+1 ; i++) {
        i_g = Nodes[i]->ind ;

        for(int j = 0 ; j < order+1 ; j++) {
            j_g = Nodes[j]->ind ;

            for(int k=0; k < quad_pts; k++) {
                quadrature(k);
                A(i_g,j_g) += w*stiffness*dshape(i,xi)*dshape(j,xi) ;
            }

            A(i_g,j_g) *= jac_det ;
        }

        for(int k=0; k < quad_pts; k++) {
                quadrature(k);
                b(i_g) += w*force*shape(i,xi) ;
        }

        b(i_g) *= jac_det ;
    }
}
