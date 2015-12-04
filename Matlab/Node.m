classdef Node
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coords %(x,y) coordinate of node
        number
        Nodes
        Elements
        dirichlet = [] ;
    end
    
    methods
        function obj = Node(n,x,y)
            obj.coords = [x; y] ;
            obj.number = n ;
        end
    end
    
end