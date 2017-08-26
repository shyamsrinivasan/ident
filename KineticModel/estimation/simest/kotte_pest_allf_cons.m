function [cons,x,p,vareps,p_usl,ac,flux] = kotte_pest_allf_cons()
x = casadi.SX.sym('x',3,1);
p = casadi.SX.sym('p',13,1);
p_usl = casadi.SX.sym('p_usl',3,1);
ac = casadi.SX.sym('ac',1,1);
flux = casadi.SX.sym(6,1);
vareps = casadi.SX.sym('vareps',2,1);

p_all = [p;p_usl;ac];

% parameters
K1ac = p_all(1);    % or 0.02
K3fdp = p_all(2);
L3fdp = p_all(3)*1e6;
K3pep = p_all(4);
K2pep = p_all(5);
vemax = p_all(6);        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = p_all(7);        % or 0.45
ne = p_all(8);             % or 2
d = p_all(9);
V4max = p_all(10);
k1cat = p_all(11);   
V3max = p_all(12);    
V2max = p_all(13);   
ac = p_all(17);

fx1 = k1cat.*x(3).*ac - flux(1).*(ac+K1ac);
fx2 = -vemax + (vemax - flux(2)).*(1+(KeFDP./x(2)).^ne);
ratio = 1+x(2)./K3fdp;
fx3 = V3max.*(ratio-1).*(ratio).^3 - flux(3).*(ratio.^4+L3fdp.*(1+x(1)./K3pep).^(-4));
fx4 = V2max.*x(1) - flux(4).*(x(1)+K2pep);
fx5 = V4max.*x(1) - flux(5);
fx6 = d.*x(3)-flux(6);

% nle
nlerhs = [fx1-fx4-fx5;...
          fx4-fx3;...
          fx2-fx6];  
      
% concentration noise cons
x_noise_lb = x - casadi.DM(xexp)*(1+vareps(1));
x_noise_ub = -x + casadi.DM(xexp)*(1-vareps(1));
      
% flux noise cons
v_noise_lb = flux - data.vexp*(1+vareps(2));
v_noise_ub = -flux + data.vexp*(1-vareps(2));      
      
cons = [fx1;fx2;fx3;fx4;fx5;fx6;nlerhs;...
        x_noise_lb;x_noise_ub;...
        v_noise_lb;v_noise_ub];





