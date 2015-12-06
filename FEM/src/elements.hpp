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
        double jac_det ;
        virtual vector<double> quadrature(int n);
        void kf_calc(MatrixXd &k_global, VectorXd &f_global);
};

class Element1d : public Element {
    private:
        Element1d(int ind_init, int order_init, int quad_init, vector<Node1d*> Nodes_init, vector<Line1d*> Edges_init);
    public:
        vector<Node1d*> Nodes;
        vector<Line1d*> Edges;
        //MatrixXd k_local;
        //VectorXd f_local;  
        virtual vector<double> quadrature(int n);
        double master_2_global( double xi ) ;
        void jacobian_calc() ;
};

#endif
