A = [-1,-0.75;1,0];
B = [1;0];
C = [1,1];
D = 0;
Plant = ss(A,B,C,D);

N = 250;
t = linspace(0,25,N);
u = [ones(N/2,1); zeros(N/2,1)];
x0 = [1;2];
[y,t,x] = lsim(Plant,u,t,x0);

figure
plot(t,y);
title('Output');
