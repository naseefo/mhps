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

load 'Cent_acc_00.txt' %Ensure the data starts from "zero".
load 'Cent_acc_90.txt'
XG = Cent_acc_00./100; % '/100' used to converted cm/s2 to m/s2
YG = Cent_acc_90./100; % '/100' used to converted cm/s2 to m/s2

NDT = length(XG);
DTN = 0.02; % Time step in seconds
NDIV = 100;   % Number of intervals the specified time interval is divided
NIT = 2;     % Number of time gauss siedel method is applied to determine iterative solution for coupled equations
NST= length(L_mass); % Number of storey
N = 1; % Number of set of runs for different parameters to be executed

for KN = 1:N
    
    %% Superstructure Iterative Paramters Intialization
    TX1 = 0.5;
    RTYTX = 1.0;
    
    %% Base Isolator Iterative Parameters Initialization
    RMBM = 1;
    TBX = 1;
    TBY = 1;
    ZETABX = 0.0;
    ZETABY = 0.0;
    BM = RMBM*L_mass(end);
    AMU = 0.1;
    
    %% Base Isolator Properties and Other Variables Calculations
    WBX = 2*pi/TBX;
    WBY = 2*pi/TBY;
    TM = BM+sum(L_mass);
    QX = 9.81*AMU*TM;
    QY = QX;
    CKABX = TM*WBX^2;;
    CKABY = TM*WBY^2;
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
    [MSB, CSB, KSB, MSBG,ISB] = MCKXY(MSBx, CSBx, KSBx, MSBGx, MSBy, CSBy, KSBy, MSBGy)
    
    % Initializations
    ID = 1;
    FAB(1:2,1) = 0.0;
    PRES(1:2,1) = 0.0;
    
    D(1:NDOFC,1:2) = 0.0;
    V(1:NDOFC,1:2) = 0.00001;
    A(1:NDOFC,1:2) = 0.0;
    
    AA(1:NDOFC,1) = 0.0;
    PD(1:NDOFC,1) = 0.0;
    PA(1:NDOFC,1) = 0.0;
    
    DD(1:NDOFC,1) = 0.0;
    DV(1:NDOFC,1) = 0.0;
    P(1:NDOFC,1:2) = 0.0;
    PINC(1:NDOFC,1) = 0.0;
    
    TIME = 0.0;
    
    for K =1:1 %NDT-1
        
        P = -1*MSBG*ISB*[XG(K) XG(K+1); YG(K) YG(K+1)];
        PINC = (P(:,2) - P(:,1))./NDIV;
        
        for KS = 1:1 %NDIV
            TIME = TIME + DT;
            P(1:NDOFC,2) = P(1:NDOFC,1) + PINC;
            EQ = -1*MSBG(1:2,1:2)*P(1:2,2);
            
            if (K+KS-2) == 0
                A(1:NSTC,2) = ACCNS(D(:,2), V(:,2), MSS, CSS, KSS, P(:,1), NSTC);
            end
            A(1:NSTC,1) = A(1:NSTC,2);  % This part can be modified by replacing A(1:NSTC,2) with A(1:NSTC,1) in the above if-statement????
            A(1:NSTC,2) = 0;
            
            ID = 0;
            if ID == 1
                EFP = EFFP(V(:,1),A(:,1),MSS,CSS,PINC,DT,NSTC);
                EFK = EFFK(MSS,CSS,KSS,DT,NSTC);                                   % NSTC input argument can be removed
                DD(1:NSTC) = inv(EFK)*EFP;                                                 % There is concern here with regards to the condition of the matrix... all the rows should have same values but there is sligh variation both in FORTRAN and MATLAB
                
                D(1:NSTC,2) = D(1:NSTC,1) + DD(1:NSTC);
                DV(1:NSTC) = 3*DD(1:NSTC)./DT - 3*V(1:NSTC,1) - 0.5*DT*A(1:NSTC,1);
                V(1:NSTC,2) = V(1:NSTC,1) + DV(1:NSTC);
                
                FAB = [CNx 0;0 CNy]*V(NSTC-1:NSTC,2) + [KNx 0;0 KNy]*D(NSTC-1:NSTC,2) - [CKABX 0;0 CKABY]*D(NDOFC-1:NDOFC,2) - [BM 0;0 BM]*EQ;
                [ID,FAB] = Stat1(FAB,[QX; QY]);
                
                
                if ID == 1
                    A(1:NSTC,2) = ACCNS(D(:,2), V(:,2), MSS, CSS, KSS, P(:,2), NSTC);
                else
                    A(1:NDOFC,2) = ACCNO(D(:,2),V(:,2),MSB,CSB,KSB,P(:,2),ISB,NDOFC,FAB);
                end
                [PA, PD] = PEAK(D(:,1), A(:,2), [EQ(1); EQ(2)], NDOFC, ISB, PA,PD);
                
                A(:,1) = A(:,2);
                V(:,1) = V(:,2);
                D(:,1) = D(:,2);
                P(:,1) = P(:,2);
                
            else
                
                DFAB(1) = 0.0;
                DFAB(2) = 0.0;
                
                for KIT = 1:NIT
                    
                    EFP = EFFP(V(:,1),A(:,1),MSB,CSB,PINC,DT,NDOFC);
                    EFK = EFFK(MSB,CSB,KSB,DT,NDOFC);                                   % NSTC input argument can be removed
                    DD(1:NDOFC) = inv(EFK)*EFP;                                                 % There is concern here with regards to the condition of the matrix... all the rows should have same values but there is sligh variation both in FORTRAN and MATLAB
                    
                    D(1:NDOFC,2) = D(1:NDOFC,1) + DD(1:NDOFC);
                    DV(1:NDOFC) = 3*DD(1:NDOFC)./DT - 3*V(1:NDOFC,1) - 0.5*DT*A(1:NDOFC,1);
                    V(1:NDOFC,2) = V(1:NDOFC,1) + DV(1:NDOFC);
                    
                    RVL = sqrt(V(NDOFC-1,2)^2 + V(NDOFC,2)^2)
                    
                    if RVL > 1.0e-4
                        DFAB(1) = QX*V(NDOFC-1,2)/RVL - FAB(1);
                        DFAB(2) = QY*V(NDOFC,2)/RVL - FAB(2);
                    end
                    
                end
                
                
                
            end
            
        end
        
    end
    
end
    ht = toc
    avg_time_per_run = ht/N
    %% Newmark Beta Parameters
