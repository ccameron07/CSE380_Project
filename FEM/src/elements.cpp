#include "elements.hpp"
#include "./eigen3/Eigen/Dense"
#include <vector>
#include <iostream>
#include <functional>


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

void Line1d::addNodes(int n, std::vector<Node1d> &nodes_g, std::vector<Node1d*> &nodes_e) {
     
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

Element1d::Element1d(int ind_init, int order_init, int quad_init, std::vector<Node1d*> Nodes_init, std::vector<Line1d*> Edges_init){
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
    std::vector<double> ws {2.0,1.0, 1.0,
                            0.8888888888888888888888889,
                            0.5555555555555555555555556,
                            0.5555555555555555555555556,
                            0.6521451548625461426269361,
                            0.6521451548625461426269361,
                            0.3478548451374538573730639,
                            0.3478548451374538573730639} ;
    
    std::vector<double> xis {0.0,
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

void Element1d::AbCalc(MatrixXd &A, VectorXd &b) {

    int size = A.rows() ;
    int i_g, j_g ;
    double f ;

    for(int i = 0 ; i < order+1 ; i++) {
        i_g = Nodes[i]->ind ;

        for(int j = 0 ; j < order+1 ; j++) {
            j_g = Nodes[j]->ind ;

            for(int k=0; k < quad_pts; k++) {
                quadrature(k) ;
                A(i_g,j_g) += w * stiffness(master_2_global(xi)) * dshape(i,xi) * dshape(j,xi) * jac_det ;
            }

        }

        for(int k=0; k < quad_pts; k++) {
                quadrature(k);
                b(i_g) += w * forcing(master_2_global(xi)) * shape(i,xi) * jac_det ;
        }

    }
}

Domain1d::Domain1d(int order_init, int nx_init, int quad_pts_init, double Xmin_init, double Xmax_init, double dirichlet_init) {
    order = order_init ;
    nx = nx_init ;
    quad_pts = quad_pts_init ;
    Xmin = Xmin_init ;
    Xmax = Xmax_init ;
    dirichlet = dirichlet_init ;
    Elements.reserve(nx) ;
    Nodes.reserve(order*nx+1) ;
    Edges.reserve(nx) ;
}

void Domain1d::build_elements() {
    int n ;
    double x ;
    double dx = (Xmax-Xmin)/nx ;
    std::vector<Node1d*> last_two_nodes ;
    std::vector<Line1d*> last_edge ;

    for(int i = 0 ; i < nx; i++) {
        n = Nodes.size();
        if (i == 0) {
            x = Xmin ;
            Nodes.emplace_back(Node1d(n, x, 0, 1, dirichlet)) ;
            n++ ;
        }
        x += dx ; 
        if (i == nx-1) {
            Nodes.emplace_back(Node1d(n, x, 0, 1, dirichlet)) ;
        } else {
            Nodes.emplace_back(Node1d(n, x)) ;
        }
        last_two_nodes.clear() ;
        last_two_nodes.emplace_back( &Nodes[n-1] ) ;
        last_two_nodes.emplace_back( &Nodes[n] ) ;

        Edges.emplace_back( Line1d(last_two_nodes[0] , last_two_nodes[1])) ;

        last_edge.clear() ;
        last_edge.emplace_back( &Edges[i] ) ;

        Elements.emplace_back( Element1d(i, order, quad_pts, last_two_nodes, last_edge)) ;
    }
}


void Domain1d::build_Ab(MatrixXd* A, VectorXd* b) {
    for(int i = 0; i < nx; i++) {
        Elements[i].forcing = forcing ;
        Elements[i].stiffness = stiffness ;
        Elements[i].AbCalc((*A), (*b)) ;
    }
}

void Domain1d::add_constraints(MatrixXd& A, VectorXd& b) {
    int n_nodes = Nodes.size() ;
    int ind ;
    for(int i = 0; i < n_nodes; i++){
        if( Nodes[i].boundary == 1){
            b -= A.col(Nodes[i].ind)*Nodes[i].BC ;
            b(Nodes[i].ind) = Nodes[i].BC ;
            A.col(Nodes[i].ind) = VectorXd::Zero(n_nodes) ;
            A.row(Nodes[i].ind) = VectorXd::Zero(n_nodes) ;
            A(Nodes[i].ind, Nodes[i].ind) = 1.0 ;
        }
    }
}
