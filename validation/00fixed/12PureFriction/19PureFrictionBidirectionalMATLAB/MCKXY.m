function [M, C, K, MG, I] = MCKXY(Mx, Cx, Kx, MGx, My, Cy, Ky, MGy);

[m n] = size(Mx);
for i = 1:m
    for j = 1:n
        M(2*i-1,2*j-1) = Mx(i,j);
        M(2*i,2*j) = My(i,j);
        
        C(2*i-1,2*j-1) = Cx(i,j);
        C(2*i,2*j) = Cy(i,j);
        
        K(2*i-1,2*j-1) = Kx(i,j);
        K(2*i,2*j) = Ky(i,j);
        
        MG(2*i-1,2*j-1) = MGx(i,j);
        MG(2*i,2*j) = MGy(i,j);
    end
    I(2*i-1,1) = 1;
    I(2*i-1,2) = 0;
    I(2*i,1) = 0;
    I(2*i,2) = 1;
end



end

