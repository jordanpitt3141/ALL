% Plot the convergance rates for the first-order in space central scheme
% with Euler time stepping for the Soliton solution using Jordan model

clear all;
close all;
format long;

% load results
A = load('D:\f77l\Boussq\Re__Data\o1normEti.dat');
B = load('D:\f77l\Boussq\Re__Data\o1normEt.dat');
C = load('D:\f77l\Boussq\Re__Data\o1normh.dat');
D = load('D:\f77l\Boussq\Re__Data\o1normu.dat');

n = max(size(A));

dx = A(:,1);
H_initial = A(:,2);
H_solution = B(:,2);
h = C(:,2);
u = D(:,2);

figure(1)
central = h(1:n);
x = dx(1:n);
loglog(x,central,'^r');
hold on
loglog([.002 .02 .02 .002],[0.0005 0.0005 0.005 0.0005]);
text(0.006,0.0003,'1');
text(0.025,0.0015,'1');
central = u(1:n);
loglog(x,central,'sr');
axis([10^-4 10 10^-5 10]);

% str1(1) = {'Error norms using the first-order central scheme dambreak problem'};
% str1(2) = {'water depth and velocity t = 20s'};
% title(str1)
xlabel('$Log_{10}\Delta x$');
ylabel('$Log_{10}L_1$');
%LEGEND('Water Depth','First-order','Velocity',4)

figure(2)
central = H_initial(1:n);
loglog(x,central,'^r');
hold on
loglog([.001 .01 .01 .001],[0.0001 0.0001 0.001 0.0001]);
text(0.012,0.00027,'1');
text(0.003,0.00007,'1');
central = H_solution(1:n);
loglog(x,central,'sr');
%axis([10^-4 10 10^-5 10]);
xlabel('$Log_{10}\Delta x?????$');
ylabel('$Log_{10}L_1$');

set (0,'defaulttextinterpreter','none')

laprint
