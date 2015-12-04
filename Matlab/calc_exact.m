function [exact] = calc_exact(x,y,n)
exact = zeros(size(x)) ;
for ii = 1:2:2*n
    for jj = 1:2*n
        Umn = (-4*(1-(-1)^ii)*(1-(-1)^jj))/((ii^2+jj^2)*ii*jj*pi^2) ;
        exact = exact+Umn*sin(ii*x).*sin(jj*y) ;
    end
end
