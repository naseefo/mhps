function [Fs2] = FS2(a0,b0,Total_Mass,Mu,xb,vb,ab)
global xb_tab theta_tab

e = sqrt((a0^2-b0^2)/a0^2);

theta_d0(1) = interp1(xb_tab,theta_tab,xb(1)); % in radians

theta_r_d0(1) = atan((b0/a0)*tan(theta_d0(1)));

theta_r_d1(1) = (a0*b0*(sec(theta_d0(1)))^2)/(a0^2+(b0*tan(theta_d0(1)))^2);

theta_r_d2(1) = (theta_r_d1(1))^2*(a0^2-b0^2)^2*(sin(2*theta_d0(1)))/(a0*b0);

c_d0(1) = a0*sin(theta_d0(1))*cos(theta_r_d0(1)) - b0*cos(theta_d0(1))*sin(theta_r_d0(1));

p_d0(1) = a0*sin(theta_d0(1))*sin(theta_r_d0(1)) + b0*cos(theta_d0(1))*cos(theta_r_d0(1));

c_d1(1) = a0*cos(theta_d0(1))*cos(theta_r_d0(1)) + b0*sin(theta_d0(1))*sin(theta_r_d0(1)) - p_d0(1)*theta_r_d1(1);

p_d1(1) = a0*cos(theta_d0(1))*sin(theta_r_d0(1)) - b0*sin(theta_d0(1))*cos(theta_r_d0(1)) + c_d0(1)*theta_r_d1(1);

c_d2(1) = -(c_d0(1) + p_d1(1)*theta_r_d1(1) + p_d0(1)*theta_r_d2(1)) - (a0*cos(theta_d0(1))*sin(theta_r_d0(1)) - b0*sin(theta_d0(1))*cos(theta_r_d0(1)))*theta_r_d1(1);

p_d2(1) = (c_d1(1)*theta_r_d1(1) + c_d0(1)*theta_r_d2(1) - p_d0(1)) + (a0*cos(theta_d0(1))*cos(theta_r_d0(1)) + b0*sin(theta_d0(1))*sin(theta_r_d0(1)))*theta_r_d1(1);

F_d0(1) = a0*sqrt(1-e^2*(sin(theta_d0(1)))^2);

F_d1(1) = -(a0*e)^2/F_d0(1);

theta_dot1(1) = (0.5*vb(1))/(F_d0(1)-c_d1(1));

theta_dot2(1) = ((F_d1(1)-c_d2(1))*theta_dot1(1)^2 - 0.5*ab(1))/(c_d1(1)-F_d0(1));

y_r_d2(1) = theta_dot1(1)^2*p_d2(1) + p_d1(1)*theta_dot2(1);

y_b_d2(1) = 2*y_r_d2(1);

%% Constants

g = 9.81; %m/s2

%%

Fs2 = Total_Mass*(g + y_b_d2(1))*Mu*sign(vb);

end

