function [ID] = Stat0(DDXO,DDYO, FAB)

 Q1 = FAB(1);
 Q2 = FAB(2);
 WD = Q1*DDXO + Q2*DDYO;
 
 if WD >= 0
     ID = 0;
 else
     ID = 1;
 end
    
end

