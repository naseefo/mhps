function [ID,FAB] = Stat1(Q,QXY)

ratio = sqrt(sum((Q.^2)./(QXY.^2)));

if (ratio-1) < 0
    FAB = Q;
    ID = 1;
else
    FAB = Q/ratio;
    ID = 0;
end

end

