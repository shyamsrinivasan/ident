function [ode,flux,D2FX,oderhs,x,p_var,p_other,acetate] =...
        kotteCASident(nc)

x = casadi.SX.sym('x',nc,1);
p_var = casadi.SX.sym('p_all',9,1);
p_other = casadi.SX.sum('p_other',4,1);
acetate = casadi.SX.sym('acetate',1,1);

p_all = [p_var;p_other;acetate];

flux = cell(6,1);

% parameters
K1ac = p_all(1);    % or 0.02
K3fdp = p_all(2);
L3fdp = p_all(3)*1e6;
K3pep = p_all(4);
K2pep = p_all(5);
V4max = p_all(6);
k1cat = p_all(7);   
V3max = p_all(8);    
V2max = p_all(9);
vemax = p_all(10);        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = p_all(11);        % or 0.45
ne = p_all(12);             % or 2
d = p_all(13);
ac = p_all(14);

flux{1} = k1cat.*x(3).*ac./(ac+K1ac);
flux{2} = vemax.*(1-1./(1+(KeFDP./x(2)).^ne));
ratio = 1+x(2)./K3fdp;
flux{3} = V3max.*(ratio-1).*(ratio).^3./...
            (ratio.^4+L3fdp.*(1+x(1)./K3pep).^(-4));
flux{4} = V2max.*x(1)./(x(1)+K2pep);  
flux{5} = V4max.*x(1);
flux{6} = d.*x(3);

oderhs = [flux{1} - flux{4} - flux{5};...
          flux{4} - flux{3};...
          flux{2} - flux{6}];
ode = casadi.Function('ode',{x,p_var,p_other,acetate},{oderhs}); 

fluxeq = [flux{1};flux{2};flux{3};flux{4};flux{5};flux{6}];
flux = casadi.Function('flux',{x,p_var,p_other,acetate},{fluxeq});
D2FX = [];





