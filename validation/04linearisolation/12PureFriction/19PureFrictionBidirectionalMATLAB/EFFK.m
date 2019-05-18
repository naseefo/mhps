function [EFK] = EFFK(M,C,K,T,N)

EFK = K + (3/T)*C + (6/T^2)*M;

end

