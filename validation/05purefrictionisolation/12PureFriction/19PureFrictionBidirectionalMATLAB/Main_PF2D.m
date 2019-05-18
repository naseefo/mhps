%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         Program to find bi-directional response of shear type                                                                                                %
%                                                                             building base-isolated with pure-friction system                                                                                                     %
%                                                                                                                                                                                                                                                             %
%                                                                                                     Developed by                                                                                                                                  %
%                                                                              Naseef Ummer, Doctoral Research Scholar, IIT Delhi                                                                                            %
%                                                                           Dr. Vasant A. Matsagar, Associate Professor, IIT Delhi                                                                                           %
%                                                                                                                                                                                                                                                             %
%     Last Updated On: 29.01.2017                                                                                                                                                                                                        %
%     Approval Status: Pending | Validated but supervisor approval required                                                                                                                                        %
%                                                                                                                                                                                                                                                             %
%     Note: Pleasre refer for variable definitions towards the end of the code                                                                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
clc

tic

format short e
format compact

%% Constants
g = 9.8; %m/s2
rad=pi/180; % for converting degrees to radians

%% Checkboxes
G_f_k_adjust=1; %Change Line 26 accordingly (L_stiffness)
G_f_c_adjust=1; %Change Line 27, 28 accordingly (L_damping ratios and G_f_damping)

%% Fixed Parameters Intialization

fk= 1;      %Common factor for L_stiffness
fm= 1;      %Common factor for L_mass
fdr=0.05;   % Common factor for L_dampring_ratios

L_mass = fm*[1 1 1 1 1]; % Enter in unit kg
L_stiffness = fk*[1 1 1 1 1]; % Enter in unit N/m
L_damping_ratios = fdr*[1 1 1 1 1]; %Specify if classical caughey damping matrix is to be constructed

G_f_damping=[0]; %Specify explicity if required but ensure G_f_c_adjust=0

% load 'AIMP_A0p5g_0p398.txt' %Ensure the data starts from "zero".
load 'Cent_acc_0.txt'
load 'Cent_acc_90.txt'
% XG =AIMP_A0p5g_0p398./2; % '/100' used to converted cm/s2 to m/s2
XG = Cent_acc_0./100;
YG = Cent_acc_90./100; % '/100' used to converted cm/s2 to m/s2

NDT = min(length(XG),length(YG));
DTN = 0.02; % Time step in seconds
NDIV = 10;   % Number of intervals the specified time interval is divided
NIT = 3;     % Number of time gauss siedel method is applied to determine iterative solution for coupled equations
NST= length(L_mass); % Number of storey
N = 1; % Number of set of runs for different parameters to be executed

for KN = 1:N
    
    %% Superstructure Iterative Paramters Intialization
    TX1 = 0.5;
    RTYTX = 1.0;
    
    %% Base Isolator Iterative Parameters Initialization
    RMBM = 1;
    TBX = 50.0;
    TBY = 50.0;
    ZETABX = 0.0;
    ZETABY = 0.0;
    BM = RMBM*L_mass(end);
    AMU = 0.01;
    
    %% Base Isolator Properties and Other Variables Calculations
    WBX = 2*pi/TBX;
    WBY = 2*pi/TBY;
    TM = BM+sum(L_mass);
    QX = 9.81*AMU*TM;
    QY = QX;
    if TBX > 10
        CKABX = 0;
    else
        CKABX = TM*WBX^2;
    end
    
    if TBY > 10
        CKABY = 0;
    else
        CKABY = TM*WBY^2;
    end
    CDABX = 2*ZETABX*WBX*TM;
    CDABY = 2*ZETABY*WBY*TM;
    DT = DTN/NDIV;
    TY1 = RTYTX*TX1;
    
    NSTC = 2*NST;
    NDOF = NST + 1;
    NDOFC = NSTC + 2;
    
    [MSSx, CSSx, KSSx, MSGx, ISSx] = MCKSS(L_mass, L_damping_ratios, L_stiffness,G_f_k_adjust,TX1,G_f_c_adjust);
    [MSSy, CSSy, KSSy, MSGy, ISSy] = MCKSS(L_mass, L_damping_ratios, L_stiffness,G_f_k_adjust,TY1,G_f_c_adjust);
    [MSS, CSS, KSS, MSG, ISS] = MCKXY(MSSx, CSSx, KSSx, MSGx, MSSy, CSSy, KSSy, MSGy);
    
    [MSBx, CSBx, KSBx, MSBGx, CNx, KNx, ISBx] = MCKSB(MSSx, CSSx, KSSx, BM, CDABX, CKABX, NST, NDOF);
    [MSBy, CSBy, KSBy, MSBGy, CNy, KNy, ISBy] = MCKSB(MSSy, CSSy, KSSy, BM, CDABY, CKABY, NST, NDOF);
    [MSB, CSB, KSB, MSBG,ISB] = MCKXY(MSBx, CSBx, KSBx, MSBGx, MSBy, CSBy, KSBy, MSBGy);
    
    % Initializations
    ID = 1;
    FAB(1:2,1) = 0.0;
    PRES(1) = 0.0;
    PRES(2) = 0.0;
    
    % X-AXIS
    DX1(1:NDOF,1) = 0.0;
    VX1(1:NDOF,1) = 1e-10;;
    AX1(1:NDOF,1) = 0.0;
    
    DX2(1:NDOF,1) = 0.0;
    VX2(1:NDOF,1) = 1e-10;
    AX2(1:NDOF,1) = 0.0;
    
    AAX(1:NDOF,1) = 0.0;
    PDX(1:NDOF,1) = 0.0;
    PAX(1:NDOF,1) = 0.0;
    
    DDX(1:NDOF,1) = 0.0;
    DVX(1:NDOF,1) = 0.0;
    PX1(1:NDOF,1) = 0.0;
    PX2(1:NDOF,1) = 0.0;
    PINCX(1:NDOF,1) = 0.0;
    
    % Y-AXIS
    DY1(1:NDOF,1) = 0.0;
    VY1(1:NDOF,1) = 1e-10;
    AY1(1:NDOF,1) = 0.0;
    
    DY2(1:NDOF,1) = 0.0;
    VY2(1:NDOF,1) = 1e-10;
    AY2(1:NDOF,1) = 0.0;
    
    AAY(1:NDOF,1) = 0.0;
    PDY(1:NDOF,1) = 0.0;
    PAY(1:NDOF,1) = 0.0;
    
    DDY(1:NDOF,1) = 0.0;
    DVY(1:NDOF,1) = 0.0;
    PY1(1:NDOF,1) = 0.0;
    PY2(1:NDOF,1) = 0.0;
    PINCY(1:NDOF,1) = 0.0;
    
    TIME = 0.0;
    
    for K =1:NDT-1
        
        if mod(K,50) == 0
            Progess = K/NDT*100;
        end
        
        PX1 = -1*diag(MSBx)*XG(K);
        PX2 = -1*diag(MSBx)*XG(K+1);
        PY1 = -1*diag(MSBy)*YG(K);
        PY2 = -1*diag(MSBy)*YG(K+1);
        PINCX = (PX2-PX1)./NDIV;
        PINCY = (PY2-PY1)./NDIV;
        
        for KS = 1:NDIV
        
            TIME = TIME + DT;
            PX2 = PX1 + PINCX;
            PY2 = PY1 + PINCY;
            EQX = -PX2(1)./MSBGx(1,1);
            EQY = -PY2(1)./MSBGy(1,1);
        
            if (K+KS-2) == 0
                AX2(1:NST) = ACCNS(DX2,VX2,MSSx,CSSx,KSSx,PX1,NST);
                AY2(1:NST) = ACCNS(DY2,VY2,MSSy,CSSy,KSSy,PY1,NST);
            end
            AX1(1:NST) = AX2(1:NST);           % This part can be modified by replacing AX2 and AY2 with AX1 and AY1, respectively in the above if-statement????
            AX2(1:NST) = 0;
            AY1(1:NST) = AY2(1:NST);
            AY2(1:NST) = 0;
            
%              ID = 0;
            if ID == 1
                ID1 = 1;
                EFPX = EFFP(VX1,AX1,MSSx,CSSx,PINCX,DT,NST);
                EFKX = EFFK(MSSx,CSSx,KSSx,DT,NST);                                   % NST input argument can be removed
                DDX(1:NST) = inv(EFKX)*EFPX;                                                  % There is concern here with regards to the condition of the matrix... all the rows should have same values but there is sligh variation both in FORTRAN and MATLAB
                
                EFPY = EFFP(VY1,AY1,MSSy,CSSy,PINCY,DT,NST);
                EFKY = EFFK(MSSx,CSSx,KSSx,DT,NST);                                   % NST input argument can be removed
                DDY(1:NST) = inv(EFKY)*EFPY;                                                  % There is concern here with regards to the condition of the matrix... all the rows should have same values but there is sligh variation both in FORTRAN and MATLAB
                
                DX2(1:NST) = DX1(1:NST) + DDX(1:NST);
                DVX(1:NST) = 3*DDX(1:NST)./DT - 3*VX1(1:NST) - 0.5*DT*AX1(1:NST);
                VX2(1:NST) = VX1(1:NST) + DVX(1:NST);
                
                DY2(1:NST) = DY1(1:NST) + DDY(1:NST);
                DVY(1:NST) = 3*DDY(1:NST)./DT - 3*VY1(1:NST) - 0.5*DT*AY1(1:NST);
                VY2(1:NST) = VY1(1:NST) + DVY(1:NST);
                
                FAB(1) = CNx*VX2(NST) + KNx*DX2(NST) - CKABX*DX2(NDOF) - BM*EQX;
                FAB(2) = CNy*VY2(NST) + KNy*DY2(NST) - CKABY*DY2(NDOF) - BM*EQY;
                
                [ID,FAB] = Stat1(FAB,[QX; QY]);
                
                
                if ID == 1
                    ID1 = 11;
                    AX2(1:NST) = ACCNS(DX2,VX2,MSSx,CSSx,KSSx,PX2,NST);
                    AY2(1:NST) = ACCNS(DY2,VY2,MSSy,CSSy,KSSy,PY2,NST);
                else
                    ID1 = 10;
                    AX2(1:NDOF) = ACCNO(DX2,VX2,MSBx,CSBx,KSBx,PX2,ISBx, NDOF,FAB(1));
                    AY2(1:NDOF) = ACCNO(DY2,VY2,MSBy,CSBy,KSBy,PY2,ISBy, NDOF,FAB(2));
                    
%                     COSTAX2NDOF = acos(AX2(NDOF)/sqrt(AX2(NDOF)^2 + AY2(NDOF)^2))*180/pi
%                     COSTFAB = acos(FAB(1)/sqrt(FAB(1)^2 + FAB(2)^2))*180/pi
                end
                    [PAX, PDX, AAX] = PEAK(DX1, AX2, EQX, NDOF, ISBx,PAX,PDX);
                    [PAY, PDY, AAY] = PEAK(DY1, AY2, EQY, NDOF, ISBy,PAY,PDY);
                    
                    AX1 = AX2;
                    VX1 = VX2;
                    DX1 = DX2;
                    PX1 = PX2;
                    
                    AY1 = AY2;
                    VY1 = VY2;
                    DY1 = DY2;
                    PY1 = PY2;
            else
                ID1 = 2;
                DFAB(1) = 0.0;
                DFAB(2) = 0.0;
                
                for KIT = 1:NIT
                    EFPX = EFFP(VX1,AX1,MSBx,CSBx,PINCX,DT,NDOF);
                    EFPX(NDOF) = EFPX(NDOF) - DFAB(1);
                    EFKX = EFFK(MSBx,CSBx,KSBx,DT,NDOF);                                   % NST input argument can be removed
                    DDX(1:NDOF) = inv(EFKX)*EFPX;                                                  % There is concern here with regards to the condition of the matrix... all the rows should have same values but there is sligh variation both in FORTRAN and MATLAB
                    
                    EFPY = EFFP(VY1,AY1,MSBy,CSBy,PINCY,DT,NDOF);
                    EFPY(NDOF) = EFPY(NDOF) - DFAB(2);
                    EFKY = EFFK(MSBx,CSBx,KSBx,DT,NDOF);                                   % NST input argument can be removed
                    DDY(1:NDOF) = inv(EFKY)*EFPY;                                                  % There is concern here with regards to the condition of the matrix... all the rows should have same values but there is sligh variation both in FORTRAN and MATLAB
                    
                    DX2(1:NDOF) = DX1(1:NDOF) + DDX(1:NDOF);
                    DVX(1:NDOF) = 3*DDX(1:NDOF)./DT - 3*VX1(1:NDOF) - 0.5*DT*AX1(1:NDOF);
                    VX2(1:NDOF) = VX1(1:NDOF) + DVX(1:NDOF);
                    
                    DY2(1:NDOF) = DY1(1:NDOF) + DDY(1:NDOF);
                    DVY(1:NDOF) = 3*DDY(1:NDOF)./DT - 3*VY1(1:NDOF) - 0.5*DT*AY1(1:NDOF);
                    VY2(1:NDOF) = VY1(1:NDOF) + DVY(1:NDOF);
                    
%                     XBR2 = sqrt(DX2(NDOF)^2 + DY2(NDOF)^2);
%                     VBR2 = sqrt(VX2(NDOF)^2 + VY2(NDOF)^2);
%                     COSTV = acos(VX2(NDOF)/VBR2)*180/pi
%                     COSTX = acos(DX2(NDOF)/XBR2)*180/pi
                    RVL = sqrt(VX2(NDOF)^2 + VY2(NDOF)^2);
                    
                    if RVL > 1.0e-5
                        DFAB(1) = QX*VX2(NDOF)/RVL - FAB(1);
                        DFAB(2) = QY*VY2(NDOF)/RVL - FAB(2);
                    end
                end
                if NIT > 1
                    FAB(1) = FAB(1) + DFAB(1);
                    FAB(2) = FAB(2) + DFAB(2);
                end
%                 Z1 = DDX(NDOF)
%                 Z2 = DDY(NDOF)
%                 Z3 = FAB
                ID = Stat0(DDX(NDOF),DDY(NDOF), FAB);
                
                if ID == 1
                    ID1 = 21;
                    VX2(NDOF) = 1e-10;
                    AX2(NDOF) = 0.0;
                    VY2(NDOF) = 1e-10;
                    AY2(NDOF) = 0.0;
                    AX2(1:NST) = ACCNS(DX2,VX2,MSSx,CSSx,KSSx,PX2,NST);
                    AY2(1:NST) = ACCNS(DY2,VY2,MSSy,CSSy,KSSy,PY2,NST);
                else
                    ID1 = 20;
                    AX2(1:NDOF) = ACCNO(DX2,VX2,MSBx,CSBx,KSBx,PX2,ISBx, NDOF,FAB(1));
                    AY2(1:NDOF) = ACCNO(DY2,VY2,MSBy,CSBy,KSBy,PY2,ISBy, NDOF,FAB(2));
                    
%                     COSTAX2NDOF = acos(AX2(NDOF)/sqrt(AX2(NDOF)^2 + AY2(NDOF)^2))*180/pi
%                     COSTFAB = acos(FAB(1)/sqrt(FAB(1)^2 + FAB(2)^2))*180/pi
                end
                [PAX, PDX, AAX] = PEAK(DX1, AX2, EQX, NDOF, ISBx,PAX,PDX);
                [PAY, PDY, AAY] = PEAK(DY1, AY2, EQY, NDOF, ISBy,PAY,PDY);
                
%                 if RVL < 1e-5
%                     VX2(NDOF) = 1e-5;
%                     VY2(NDOF) = 1e-5;
%                 end
                
                AX1 = AX2;
                VX1 = VX2;
                DX1 = DX2;
                PX1 = PX2;
                
                AY1 = AY2;
                VY1 = VY2;
                DY1 = DY2;
                PY1 = PY2;
                
                
                
            end
            
        end
        
        XBR = sqrt(DX2(NDOF)^2 + DY2(NDOF)^2);
        FBR = sqrt(FAB(1)^2 + FAB(2)^2);
        RESULT01(K,:) = [TIME ID1 AAX(1) AAY(1) DX2(NDOF) DY2(NDOF) FAB(1) FAB(2) XBR FBR];
        
        RESULT02(K,:) = [TIME ID1 AX1(1) VX1(1) DX1(1) PX1(1) AY1(1) VY1(1) DY1(1) PY1(1) AX2(1) VX2(1) DX2(1) PX2(1) AY2(1) VY2(1) DY2(1) PY2(1) FAB(1) FAB(2)];
        RESULT03(K,:) = [TIME ID1 AX1(6) VX1(6) DX1(6) PX1(6) AY1(6) VY1(6) DY1(6) PY1(6) AX2(6) VX2(6) DX2(6) PX2(6) AY2(6) VY2(6) DY2(6) PY2(6) FAB(1) FAB(2)];
        
    end
    
end

ht = toc
avg_time_per_run = ht/N

