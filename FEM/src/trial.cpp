#include "elements.hpp"
#include "./eigen3/Eigen/Dense"
#include <vector>
#include <iostream>

class dummy {
	public:
		std::vector<Node> Nodes ;
		void add_node(){Nodes.push_back(Node());}
} ;

int main () {
	dummy D ;
	D.add_node() ;
	std::cout << D.Nodes.size() << "   " << D.Nodes.capacity() << std::endl ;
    D.add_node() ;
	std::cout << D.Nodes.size() << "   " << D.Nodes.capacity() << std::endl ;
	{
		D.add_node() ;
		D.add_node() ;
		D.add_node() ;
		std::cout << D.Nodes.size() << "   " << D.Nodes.capacity() << std::endl ;
	}
	std::cout << D.Nodes.size() << "   " << D.Nodes.capacity() << std::endl ;
	return 0 ;
}