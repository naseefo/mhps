

clc
clear

Ref = 2;
RSN = 179;


x = load(['EQX-' num2str(Ref) '.csv']);
y = load(['EQY-' num2str(Ref) '.csv']);

figure(1)
subplot(3,1,1)
plot(x(:, 1), x(:, 2), y(:, 1), y(:, 2))
xlabel('Time (sec)')
ylabel('Ground acceleration (g)')
grid on
legend('X' , 'Y')

subplot(3,1,2)
plot(x(:, 1), x(:, 2))
xlabel('Time (sec)')
ylabel('Ground acceleration (g)')
grid on
legend('X')

subplot(3,1,3)
plot(y(:, 1), y(:, 2))
xlabel('Time (sec)')
ylabel('Ground acceleration (g)')
grid on
legend('Y')


copyfile(['EQX-' num2str(Ref) '.csv'], ['PEER_EQX-' num2str(RSN) '_FN.csv']);
copyfile(['EQY-' num2str(Ref) '.csv'], ['PEER_EQY-' num2str(RSN) '_FP.csv']);