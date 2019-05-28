
clc
clear

Mb=1;
g=9.8;

a0=0.2;

% b0=0.2;
% e=sqrt((a0^2-b0^2)/(a0^2))

e = 0.5;
b0 = sqrt(a0^2-(e*a0)^2);

xb=0;
Fb=0;
theta_rot=0;
val1=0;

index=0;
F = @(x)sqrt(1-(e.^2)*(sin(x)).^2);

for theta_r=0:0.01:pi/2
    index=index+1;
    
    theta=atan((a0/b0)*tan(theta_r));
   
    p=a0*sin(theta)*sin(theta_r)+b0*cos(theta)*cos(theta_r);
    c=a0*sin(theta)*cos(theta_r)-b0*cos(theta)*sin(theta_r);    
    I_theta=quad(F,0,theta);
    
    
    theta_rot(index)=theta_r;
    xb(index)=2*a0*I_theta-2*c;
    Fb(index)=Mb*g*(c/p);
    
    val1(index)=xb(index)/2+c;
    val2(index)=p;
    val3(index)=c;
    
    temp1_xb(index) = xb(index);
    temp1_theta_r(index) = theta_r;
end

xb=xb';
Fb=Fb';
theta_rot=theta_rot';

figure(1)
plot(xb,Fb)
% axis([0 0.7 0 1]);

figure(2)
plot(temp1_theta_r.*180/pi,temp1_xb)