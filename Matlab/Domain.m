classdef Domain
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X
        Y%(x,y) coordinate of node
        nx
        ny
        order
        dirichlet = 0 ;
        N
        E
        node_index
    end
    
    methods
        
        function obj = buildNodes(obj)
            n = 1 ;
            dx = (obj.X(2)-obj.X(1))/obj.nx ;
            dy = (obj.Y(2)-obj.Y(1))/obj.ny ;
            
            obj.N = Node.empty(obj.nx+obj.ny+2,0) ;
            obj.node_index = zeros(obj.nx+1, obj.ny+1) ;
            
            for ii = 1:obj.nx+1
                for jj = 1:obj.ny+1
                    obj.N(n) = Node(n, obj.X(1)+dx*(ii-1), obj.Y(1)+dy*(jj-1)) ;
                    obj.node_index(ii,jj) = n ;
                    if (ii == 1 || ii == obj.nx+1 || jj == 1 || jj == obj.ny+1)
                        obj.N(n).dirichlet = obj.dirichlet ;
                    end
                    n = n+1 ;
                end
            end
        end
        
        function obj = buildElements(obj)
            n = 1 ;
            obj.E = Element.empty ;
            for ii = 1:obj.nx
                for jj = 1:obj.ny
                    %Build lower left triangular element
                    
                    obj.E(n).Nodes(1) = obj.N(obj.node_index(ii,jj)) ;
                    obj.E(n).Nodes(2) = obj.N(obj.node_index(ii+1,jj)) ;
                    obj.E(n).Nodes(3) = obj.N(obj.node_index(ii,jj+1)) ;
                    
                    n = n+1 ;
                    %Build upper right triangular element
                    
                    obj.E(n).Nodes(1) = obj.N(obj.node_index(ii+1,jj)) ;
                    obj.E(n).Nodes(2) = obj.N(obj.node_index(ii+1,jj+1)) ;
                    obj.E(n).Nodes(3) = obj.N(obj.node_index(ii,jj+1)) ;
                    
                    n = n+1 ;
                end
            end
        end
        
        function [kg, fg] = buildkf(obj)
            kg = zeros(length(obj.N)) ;
            fg = zeros(length(obj.N),1) ;
            for ii = 1:length(obj.E)
                obj.E(ii) = obj.E(ii).calc_jacobian ;
                obj.E(ii) = obj.E(ii).buildk ;
                obj.E(ii) = obj.E(ii).buildf ;
                kg = kg + obj.E(ii).scatterk(length(obj.N)) ;
                fg = fg + obj.E(ii).scatterf(length(obj.N)) ;
            end
        end
        
        function [kg, fg] = applyBC(obj, kg, fg)
            for ii = 1:length(obj.N)
                if isempty(obj.N(ii).dirichlet)
                else
                    fg = fg+kg(:,obj.N(ii).number)*obj.N(ii).dirichlet ;
                    fg(obj.N(ii).number) = obj.N(ii).dirichlet ;
                    kg(:,obj.N(ii).number) = 0 ;
                    kg(obj.N(ii).number,:) = 0 ;
                    kg(obj.N(ii).number,obj.N(ii).number) = 1 ;
                end
            end
        end
        
        function domainPlot(obj)
            figure
            for ii = 1:length(obj.E)
                hold on
                obj.E(ii).endpt_plot(ii) ;
            end
        end
        
        function [x, y] = getCoords(obj)
            x = zeros(size(obj.N)) ;
            y = zeros(size(obj.N)) ;
            for ii = 1:length(obj.N)
                tmp = obj.N(ii).coords() ;
                x(ii) = tmp(1) ;
                y(ii) = tmp(2) ;
            end
        end
        

    end
    
end