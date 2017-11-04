function dx = g1sNLAE(x,model,p)
% x1 = pRB
% x2 = e2f1;
% x3 = cycdi
% x4 = cycda;
% x5 = ap1
% x6 = pRBp
% x7 = pRBpp
% x8 = cycEi
% x9 = cycEa

% parameters - Swat, et al., 2004, Bioinformatics
if nargin<3
    k1 = 1;k2 = 1.6;k3 = .05;
    k16 = .4;k61 = .3;k34 = .04;k43 = .01;k67 = .7;k76 = .1;
    k23 = .3;k25 = .9;k28 = .06;
    k89 = .07;k98 = .01;
    a = .04;
    j11 = .5;j12 = 5;j13 = .002;j15=.001;j18 = .6;
    j61 = 5;j62 = 8;j63 = 2;j65 = 6;j68 = 7;
    Km1 = .5;Km2 = 4;Km4 = .3;Km9 = .005;kp = .05;
    phipRB = .005;phiE2F1 = .1;phiCycDi = .023;phiCycDa = .03;
    phiAP1 = .01;phipRBp = .06;phipRBpp = .04;phiCycEi = .06;phiCycEa = .05;
    Fm = .006;
else
    k1 = p(1);k2 = p(2);k3 = p(3);
    k16 = p(4);k61 = p(5);k34 = p(6);k43 = p(7);k67 = p(8);k76 = p(9);
    k23 = p(10);k25 = p(11);k28 = p(12);
    k89 = p(13);k98 = p(14);
    a = p(15);
    j11 = p(16);j12 = p(17);j13 = p(18);j15 = p(19);j18 = p(20);
    j61 = p(21);j62 = p(22);j63 = p(23);j65 = p(24);j68 = p(25);
    Km1 = p(26);Km2 = p(27);Km4 = p(28);Km9 = p(29);kp = p(30);
    phipRB = p(31);phiE2F1 = p(32);phiCycDi = p(33);phiCycDa = p(34);
    phiAP1 = p(35);phipRBp = p(36);phipRBpp = p(37);phiCycEi = p(38);
    phiCycEa = p(39);Fm = p(40);
end

dx = zeros(6,1);
dx = cons(dx,x);
% core module only
% dx(1) = k1*(x(2)/(Km1+x(2)))*j11/(j11+x(1)) - phipRB*x(1);
% 
% dx(2) = kp + (k2*(a^2+x(2)^2)/(Km2^2+x(2)^2))*j12/(j12+x(1)) -...
%         phiE2F1*x(2);

% core + extension 1
dx(1) = k1*x(2)/(Km1+x(2))*j11/(j11+x(1))*j61/(j61+x(6)) -...
        k16*x(1)*x(4) + k61*x(6) - phipRB*x(1);

dx(2) = kp + k2*(a^2+x(2)^2)/(Km2^2+x(2)^2)*j12/(j12+x(1))*j62/(j62+x(6)) -...
        phiE2F1*x(2);
    
dx(3) = k3*x(5) + k23*x(2)*j13/(j13+x(1))*j63/(j63+x(6)) + k43*x(4) -...
        k34*x(3)*x(4)/(Km4+x(4)) - phiCycDi*x(3);
    
dx(4) = k34*x(3)*x(4)/(Km4+x(4)) - k43*x(4) - phiCycDa*x(4);

dx(5) = Fm + k25*x(2)*j15/(j15+x(1))*j65/(j65+x(6)) - phiAP1*x(5);
dx(6) = k16*x(1)*x(4) - k61*x(6) - phipRBp*x(6);

% dx(1) = k1*x(2)/(Km1+x(2))*j11/(j11+x(1))*j61/(j61+x(6)) -...
%         k16*x(1)*x(4) + k61*x(6) - phipRB*x(1);
% 
% dx(2) = kp + k2*(a^2+x(2)^2)/(Km2^2+x(2)^2)*j12/(j12+x(1))*j62/(j62+x(6)) -...
%         phiE2F1*x(2);
%     
% dx(3) = k3*x(5) + k23*x(2)*j13/(j13+x(1))*j63/(j63+x(6)) + k43*x(4) -...
%         k34*x(3)*x(4)/(Km4+x(4)) - phiCycDi*x(3);
%     
% dx(4) = k34*x(3)*x(4)/(Km4+x(4)) - k43*x(4) - phiCycDa*x(4);
% 
% dx(5) = Fm + k25*x(2)*j15/(j15+x(1))*j65/(j65+x(6)) - phiAP1*x(5);
% dx(6) = k16*x(1)*x(4) - k61*x(6) - k67*x(6)*x(9) + k76*x(7) - phipRBp*x(6);
% dx(7) = k67*x(6)*x(9) - k76*x(7) - phipRBpp*x(7);
% dx(8) = k28*x(2)*j18/(j18+x(1))*j68/(j68+x(6)) +...
%         k98*x(9) - k89*x(8)*x(9)/(Km9+x(9)) - phiCycEi*x(8);
% dx(9) = k89*x(8)*x(9)/(Km9+x(9)) - k98*x(9) - phiCycEa*x(9);