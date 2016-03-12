% o2c = 250e-3; %mM
% o2e = 273e-3; %mM
% ko2e = 3.4e-6;
% ko2c = 3.4e-6;
% vo2t = (469.*o2e/ko2e)./(1+o2e/ko2e+o2c/ko2c);
% e = 23.36./vo2t;
% 
% o2e = 0:0.0001:1000e-3; %mM
% vo2t = e.*(469.*o2e/ko2e)./(1+o2e/ko2e+o2c/ko2c);
% 
% hold on
% plot(o2e,vo2t,'r-');
% plot(o2e,e,'k-');

Vmax = 2200;
K1_K2 = 900;
p3 = 10^(-2.18);
hout = 1e-5:1e-6:1e-3;
hin = 1e-4;

vsyn = Vmax*(hin./hout)./(hin./hout + K1_K2 + p3./hout);
plot(hout,vsyn);

hold on

Vmax = 2200;
Km = 6.6e-6;
Ki = 7.3e-9;
hout = 1e-5:1e-6:1e-3;
hin = 1e-4;

vsyn = Vmax*hin./(hin + Km*(1+hout./Ki));
plot(hout,vsyn);

hin = 1e-4;
hout = 1e-5:1e-6:1e-3;
atp
adp
x = hin/hout*1/K;
D = 1 + x + x^2 + x^3 + x^4;
kfwd = k2*1/(1+Kc)*1/D;
kbkw = k1*Kc/(1+Kc)*(x^4/D);

vsyn = kfwd*atp - kbkw*adp*pi;