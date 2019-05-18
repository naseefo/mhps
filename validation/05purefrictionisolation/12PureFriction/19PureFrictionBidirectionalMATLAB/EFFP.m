function [EFP] = EFFP(V1,A1,SM,CD,PINC,T,N)

EFP = PINC(1:N,1);
X1 = SM*((6/T)*V1(1:N) + 3*A1(1:N));
X2 = CD*(3*V1(1:N) + 0.5*T*A1(1:N));
EFP = EFP + X1 + X2;

end

