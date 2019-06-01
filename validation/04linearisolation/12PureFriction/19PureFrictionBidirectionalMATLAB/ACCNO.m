function [A] = ACCNO(D,V,SM,CD,SK,P,I,N,FAB)

FC = CD*V(1:N);
FD = SK*D(1:N);
EP1 = P(1:N,1) - FC - FD;

[m n] = size(FAB);

if m == 2
    EP1(N-1:N,1) = EP1(N-1:N,1) - FAB;
    A(N-1:N,1) = EP1(N-1:N,1)./diag(SM(N-1:N, N-1:N));
    A(1:N-2) = (EP1(1:N-2,1)-SM(1:N-2,1:N-2)*I(1:N-2,:)*A(N-1:N,1))./diag(SM(1:N-2,1:N-2));
else
    EP1(N) = EP1(N) - FAB;
    A(N) = EP1(N)/SM(N,N);
    A(1:N-1) = (EP1(1:N-1,1)-SM(1:N-1,1:N-1)*I(1:N-1,1)*A(N))./diag(SM(1:N-1,1:N-1));
end

end

