%
s = tf('s');
sys = 1/(s^5+5*s^4+10*s^3+10*s^2+5*s+3)
h = pzplot(sys);
grid on
