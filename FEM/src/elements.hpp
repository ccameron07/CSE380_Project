#ifndef ELEMENTS_H
#define ELEMENTS_H
#include <vector>
#include <functional>
#include "./eigen3/Eigen/Dense"

using namespace Eigen ;
using namespace std ;

class Node {
    private:
        Node(){ } 
    public:
        int ind, df ;
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

class Node2d : public Node {
    public:
        Vector2d coords ;
        //Constructor Method
        Node2d(int n_init, double x, double y, int df_init = 0, bool boundary_init = false, double BC_init = 0.0) ;
};

class Line1d {
    public:
        vector<Node1d*> nodes;
        
        Line1d(Node1d* n1, Node1d* n2) ;
        void addNodes(int n, vector<Node1d>& nodes_g, vector<Node1d*>& nodes_e) ;
};

class Element {
    protected:
        Element(){ }
    public:
        int ind, order, quad_pts ;
        double jac_det, w, xi ;
        virtual void quadrature(int n) { };
        void kfCalc(MatrixXd &A, VectorXd &b);
};

class Element1d : public Element {
    public:
        Element1d(int ind_init, int order_init, int quad_init, vector<Node1d*> Nodes_init, vector<Line1d*> Edges_init);
        vector<Node1d*> Nodes;
        vector<Line1d*> Edges;
        //MatrixXd k_local;
        //VectorXd f_local;  
        virtual void quadrature(int n);
        double master_2_global( double xi ) ;
        void jacobian_calc() ;
        void kfCalc(MatrixXd &A, VectorXd &b);
        double shape(int n, double xi);
        double dshape(int n, double xi);

};

class Domain1d {
    private:
        Domain1d() { }
    public:
        Domain1d(int order_init, int nx_init, int quad_pts_init, double Xmin_init, double Xmax_init, double dirichlet_init) ;
        int order, nx, quad_pts;
        double Xmin, Xmax, dirichlet;
        std::function<double(double)> stiffness ;
        std::function<double(double)> forcing ;
        vector<Element1d> Elements ;
        vector<Line1d> Edges ;
        vector<Node1d> Nodes ;
        void build_elements() ;
        void add_constraints() ;
} ;

#endif
