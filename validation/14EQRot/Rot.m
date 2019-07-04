

clc
clear

eq = load('EQ-4.csv');

t = eq(:,1);

u1 = eq(:,2);
alpha1 = 140;

u2 = eq(:,3);
alpha2 = 230;

strike = 10;

thetax = 0;
thetay = 90;

idx = 1;
for i = 1:1:360
    
    beta1 = (i - alpha1 - thetax)*pi/180;
    beta2 = (i - alpha2 - thetay)*pi/180;
    
    
%     u_fp = u1*cos(beta1) + u2*cos(beta2);
%     u_fn = u1*sin(beta1) + u2*sin(beta2);
    
    u_fp = u1*cos(beta1) - u2*sin(beta2);
    u_fn = u1*sin(beta1) + u2*cos(beta2);
    
    maxufp(idx) = max(abs(u_fp));
    maxufn(idx) = max(abs(u_fn));
    idx = idx + 1;
    
    sprintf('FP: %8.4f g | Paper = %8.4f g || FN: %8.4f g | Paper = %8.4f g', [max(abs(u_fp)) 0.3 max(abs(u_fn)) 0.5])
    
end

i = 0:1:360;
figure(3)
subplot(2,1,1)
plot(i, maxufp, i, maxufn)
xlabel('Strike')
ylabel('PGA')
legend('FP', 'FN')

subplot(2,1,2)
plot(maxufp, maxufn)
xlabel('PGA-FP')
ylabel('PGA-FN')


