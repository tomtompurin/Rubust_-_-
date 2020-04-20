function F=notch_f(omega,zeta_n,zeta_d)
num=[1,2*zeta_n*omega,omega^2];
den=[1,2*zeta_d*omega,omega^2];
F=ss(tf(num,den));
end