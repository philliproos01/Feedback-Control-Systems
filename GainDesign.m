A = [2 1 0; 1 0 0; 3 0 1];
B = [1; 0; 0]; 
C = [0 0 1];
D = 0;
G = ss(A,B,C,D, 0.1);
x0 = [2 ; 0; -3];  % initial state

%uncomment the poles you DONT want to use
%initial(G,x0)
%poles = [0, 0, 0];
%poles = [0.5, (0.5+0.01i), (0.5+0.01i)]
poles = [0, 0.9, -0.9];

K = acker(A, B, poles);
N = (A-B*K)

t = 0:0.1:20;
sys = ss(N,B,C,D, 0.1);
initial(sys,x0,t)

