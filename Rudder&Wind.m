%rudder angle, δ
s = tf('s');
sys_a = (s + 0.0068)/(s*(s +  0.2647)*(s + 0.0063));
K=0:0.01:100;
figure(1)
rlocus(sys_a, K)



%[k,poles] = rlocfind(sys)

%wind changes, w
s = tf('s');
sys_b = (1)/(s*(s +  0.2647)*(s + 0.0063));
K=0:0.01:100;
figure(2)
rlocus(sys_b, K)



%[k,poles] = rlocfind(sys)
