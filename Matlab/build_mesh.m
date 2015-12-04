% clear
% E(1) = Element;
% 
% E(1).Nodes(1) = Node(1, 0.0, 0.0) ;
% E(1).Nodes(2) = Node(2, 0.5, 0.0) ;
% E(1).Nodes(3) = Node(3, 0.0, 0.5) ;
% 
% E(2) = Element;
% 
% E(2).Nodes(1) = Node(2, 0.5, 0.0) ;
% E(2).Nodes(2) = Node(4, 0.5, 0.5) ;
% E(2).Nodes(3) = Node(3, 0.0, 0.5) ;
% 
% for ii = 1:length(E)
%     
%     E(ii) = E(ii).calc_jacobian ;
%     
%     det(E(ii).jacobian)
%     
% %     [~, del_psi] = E(ii).shape(3, [0.1; 0.0])
% %     [~, del_psi] = E(ii).shape(2, [0.1; 0.0])
%     
% %     [w, master] = E(ii).quad(3, 1)
%     
%     E(ii) = E(ii).buildk ;
%     E(ii) = E(ii).buildf ;
%     
%     kgl(:,:,ii) = E(ii).scatterk(4) ;
% end
% 
% %%
% clear
% nx = 4 ;
% ny = 4 ;
% 
% n = 1 ;
% Nodes = Node.empty(nx+ny+2,0) ;
% node_index = zeros(nx+1, ny+1) ;
% for ii = 1:nx+1
%     for jj = 1:ny+1
%         Nodes(n) = Node(n, ii-1, jj-1) ;
%         node_index(ii,jj) = n ;
%         n = n+1 ;
%     end
% end
% 
% n = 1 ;
% E = Element.empty ;
% for ii = 1:nx
%     for jj = 1:ny
%         %Build lower left triangular element
%         
%         E(n).Nodes(1) = Nodes(node_index(ii,jj)) ;
%         E(n).Nodes(2) = Nodes(node_index(ii+1,jj)) ;
%         E(n).Nodes(3) = Nodes(node_index(ii,jj+1)) ;
%         
%         n = n+1 ;
%         %Build upper right triangular element
%         
%         E(n).Nodes(1) = Nodes(node_index(ii+1,jj)) ;
%         E(n).Nodes(2) = Nodes(node_index(ii+1,jj+1)) ;
%         E(n).Nodes(3) = Nodes(node_index(ii,jj+1)) ;
%         
%         n = n+1 ;
%     end
% end
% 
% %%
% figure
% for ii = 1:length(E)
%     hold on
%     E(ii).endpt_plot(ii) ;
% end
% 
% %%
% kg = zeros(length(Nodes)) ;
% fg = zeros(length(Nodes),1) ;
% for ii = 1:length(E)
%     E(ii) = E(ii).calc_jacobian ;
%     E(ii) = E(ii).buildk ;
%     E(ii) = E(ii).buildf ;
%     kg = kg + E(ii).scatterk(length(Nodes)) ;
%     fg = fg + E(ii).scatterf(length(Nodes)) ;
% end


%%
D = Domain ;
D.X = [0 pi] ;
D.Y = [0 pi] ;
D.nx = 32 ;
D.ny = D.nx ;

D = D.buildNodes ;
D = D.buildElements ;
%D.domainPlot

[kg, fg] = D.buildkf ;
[kg1, fg1] = D.applyBC(kg, fg) ;
[x, y] = D.getCoords ;
z = kg1\fg1 ;
x = reshape(x,D.ny+1,D.nx+1) ;
y = reshape(y,D.ny+1,D.nx+1) ;
z = reshape(z,D.ny+1,D.nx+1) ;
%figure
%contourf(x, y, z)
exact = calc_exact(x,y,15) ;
n = [n D.nx*D.ny] ;
e = [e trapz(x(1,:),trapz(y(:,1),(exact+z).^2,2),1)] ;
    

% %%
% for ii = 1:5
%     exact(:,:,ii) = calc_exact(x,y,ii+2) ;
% end
