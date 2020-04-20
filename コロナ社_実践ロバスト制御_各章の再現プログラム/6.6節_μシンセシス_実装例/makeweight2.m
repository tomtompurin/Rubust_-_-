%% 2次の伝達関数を返す関数
% gL,gH(低域，高域のゲイン),wx(零クロス周波数)zeta(折れてん周波数の減衰比)
function W=makeweight2(gL,wx,gH,zeta)
wH=wx*sqrt(gH);
wL=wx*sqrt(gL);
num=[1,2*zeta*wL,wL^2]*gH;
den=[1,2*zeta*wH,wH^2];
W=ss(tf(num,den));
end