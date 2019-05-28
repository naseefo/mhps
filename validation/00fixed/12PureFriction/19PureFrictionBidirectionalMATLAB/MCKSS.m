function [M,C,K, MSG, I] = MCKSS(L_mass, L_damping_ratios, L_stiffness,G_f_k_adjust,TS,G_f_c_adjust)

NST = length(L_mass);

%% Assembling global fixed-base mass matrix
for i = 1:NST
    M(i,i) = L_mass(i);
end

%% Assembling global fixed-base stiffness matrix
if (NST==1)
    K(1,1)=L_stiffness(1);
else
    
    for i = NST:-1:1
        if i==1
            K(i,1)=L_stiffness(i);
            K(i,2)=(-1)*L_stiffness(i);
        end
        if (i>1 && i<NST)
            for j=(i-1):(i+1)
                if j==(i-1)
                    K(i,j)=(-1)*L_stiffness(j);
                end
                if j==i
                    K(i,j)=L_stiffness(i-1)+L_stiffness(i);
                end
                if j==(i+1)
                    K(i,j)=(-1)*L_stiffness(j);
                end
            end
        end
        
        if i==NST
            K(i,i-1)=(-1)*L_stiffness(i);
            K(i,i)=L_stiffness(i-1)+L_stiffness(i);
        end
    end
end
[f_mode_shape, f_ang_freq_sq]=eig(K,M);

%% Adjusting stiffness matrix for fixed-base fundamental time period if

if (G_f_k_adjust == 1)
    WS=2*pi/TS;
    RT=WS/(f_ang_freq_sq(1,1).^0.5);
    K= K*(RT^2);
end


%% Calculation of undamped mode shape matrix and natural time period
[f_mode_shape, f_ang_freq_sq]=eig(K,M);
f_timeperiod_superstructure=2*pi*diag(f_ang_freq_sq).^(-0.5);
f_ang_freq_superstructure=2*pi./f_timeperiod_superstructure;

%% Calculation of fixed-base global damping matrix (Caughey damping matrix)
% It is a classical damping matrix
if (G_f_c_adjust == 1)
    for i=1:NST
        for j=1:NST
            c_w(i,j)=f_ang_freq_superstructure(i)^(2*(j-1)-1);
        end
        c_cr(i,1)=2*L_damping_ratios(i);
    end
    
    c_a=inv(c_w)*c_cr;
    C=0;
    
    for i=1:NST
        C=C+c_a(i)*M*(inv(M)*K)^(i-1);
    end
end

MSG = M;
I = ones(NST,1);
% Othomass=f_mode_shape'*M*f_mode_shape;
% Othostiffness=f_mode_shape'*K*f_mode_shape;
% Othodamping=(f_mode_shape'*C*f_mode_shape);


% t1=f_ang_freq_sq.^0.5;
% ((f_mode_shape'*C*f_mode_shape)./(2*t1))./M; % SOME PROBLEM HERE IN ORTHOGONALIZING THE MASS MATRIX... Why inf is coming? Diagonal part seems correct...



end

