function [PA, PD, AA] = PEAK(D1, A2, EQ, N, I, PAT, PDT)
global PA PD
[m n] = size(EQ);


if m ==2
AA(N-1:N,1) = A2(N-1:N,1) + EQ;
AA(1:N-2,1) = A2(1:N-2,1) + I(1:N-2,:)*A2(N-1:N,1) + I(1:N-2,:)*EQ;

for i = 1:N
     if (abs(PAT(i,1)) - abs(AA(i,1))) < 0
        PA(i,1) = AA(i,1);
    else
        PA(i,1) = PAT(i,1);
    end
        
    if (abs(PDT(i,1)) - abs(D1(i,1))) < 0
        PD(i,1) = D1(i,1);
    else
        PD(i,1) = PDT(i,1);
    end 
    
end

else
AA(N,1) = A2(N,1) + EQ;
AA(1:N-1,1) = A2(1:N-1,1) + I(1:N-1,1)*A2(N,1) + I(1:N-1,1)*EQ;

for i = 1:N
    if (abs(PAT(i,1)) - abs(AA(i,1))) < 0
        PA(i,1) = AA(i,1);
    else
        PA(i,1) = PAT(i,1);
    end
        
    if (abs(PDT(i,1)) - abs(D1(i,1))) < 0
        PD(i,1) = D1(i,1);
    else
        PD(i,1) = PDT(i,1);
    end
end

end






end

