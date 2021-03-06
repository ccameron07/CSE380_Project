classdef Element
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Nodes = Node.empty %Vector of nodes associated with element
        q_pts = 3 ;
        order = 1 ;
        k = @(x) 1 ;
        %f = @(x) (sin(pi*x(1))*sin(pi*x(2))).^3 ;
        f = @(x) 1 ;
        jacobian = [];
        kl = [] ;
        fl = [] ;
    end
    
    methods
        
        function obj = Element(k, f)
            if nargin == 2
                obj.k = k ;
                obj.f = f ;
            end
        end
        
        function coords = get_coords(obj,i)
            coords = obj.Nodes(i).coords ;
        end
        
        function obj = calc_jacobian(obj)
            endpts = [obj.get_coords(1), obj.get_coords(2), obj.get_coords(3)];
            obj.jacobian = [(endpts(:,2)-endpts(:,1)), (endpts(:,3)-endpts(:,1))];
        end
            
        function local = Te(obj, master)
            if isempty(obj.jacobian)
                obj = obj.calc_jacobian ;
            end
            local = obj.jacobian*master + obj.get_coords(1) ;
        end
        
        function endpt_plot(obj,n)
            endpts = [obj.get_coords(1), obj.get_coords(2), obj.get_coords(3)];
            if mod(n,2)
                plot([endpts(1,:),endpts(1,1)],[endpts(2,:),endpts(2,1)],...
                ':o','linewidth',2)
            else
                plot([endpts(1,:),endpts(1,1)],[endpts(2,:),endpts(2,1)],...
                'r--o','linewidth',2)
            end
            
            cent = sum(endpts,2)./3 ;
                
                text(cent(1),cent(2),num2str(n))  
        end
        
        function [psi, del_psi] = shape(obj, i, master)
            switch i
                case 1
                    psi = 1-master(1)-master(2) ;
                    del_psi = [-1; -1];
                case 2
                    psi = master(1) ;
                    del_psi = [1; 0];
                case 3
                    psi = master(2) ;
                    del_psi = [0; 1];
            end
            del_psi = inv(obj.jacobian)'*del_psi ;
        end
        
        function [w, master] = quad(obj, q_pts, i)
            switch q_pts
                case 1
                    w = 0.5 ;
                    master = [1/3; 1/3] ;
                case 3
                    w = ones(1,3)*1/6 ;
                    master = [[0; 0.5] [0.5; 0], [0.5; 0.5]] ;
                    w = w(i) ;
                    master = master(:,i) ;
            end
        end
        
        
        function obj = buildk(obj)
            obj.kl = zeros(2*obj.order+1) ;
            for ii = 1:2*obj.order+1
                for jj = 1:2*obj.order+1
                    for kk = 1:obj.q_pts
                        [w, master] = obj.quad(obj.q_pts,kk);
                        [~, del_psi2] = obj.shape(jj, master);
                        [~, del_psi1] = obj.shape(ii, master);
                        obj.kl(ii,jj) = obj.kl(ii,jj)+w*obj.k(obj.Te(master))*...
                            (del_psi1'*del_psi2)*det(obj.jacobian);
                    end
                end
            end
        end
        
        function obj = buildf(obj)
            obj.fl = zeros(2*obj.order+1,1) ;
            for ii = 1:2*obj.order+1
                for kk = 1:obj.q_pts
                    [w, master] = obj.quad(obj.q_pts,kk);
                    [psi, ~] = obj.shape(ii, master);
                    obj.fl(ii) = obj.fl(ii)+w*obj.f(obj.Te(master))*...
                        psi*det(obj.jacobian) ;
                end
            end
        end

            
        function [kg] = scatterk(obj,n_nodes)
            kg = zeros(n_nodes) ;
            for ii = 1:length(obj.Nodes)
                for jj = 1:length(obj.Nodes)
                    kg(obj.Nodes(ii).number,obj.Nodes(jj).number) = obj.kl(ii,jj) ;
                end
            end
        end
        
        function [fg] = scatterf(obj,n_nodes)
            fg = zeros(n_nodes,1) ;
            for ii = 1:length(obj.Nodes)
                    fg(obj.Nodes(ii).number) = obj.fl(ii) ;
            end
        end
        
    end
end
    

