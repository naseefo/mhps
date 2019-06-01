function [A2] = ACCNS(D,V,SM,CD,SK,P,N)

FC = CD*V(1:N);
FD = SK*D(1:N);
EP1 = P(1:N,1) - FC - FD;
A2 = EP1./diag(SM);
end

