function [ode,dfx_sym,D2FX,oderhs,x,p_all,ident_c,p_useless,acetate] =...
        kotteCASident(idx,nc)

x = casadi.SX.sym('x',nc,1);
p_all = casadi.SX.sym('p_all',12,1);
ident_c = casadi.SX.sym('ident_c',1,1);
p_useless = casadi.SX.sym('p_useless',3,1);
acetate = casadi.SX.sym('acetate',1,1);

p = [p_all(1:idx-1);ident_c;p_all(idx:end);p_useless;acetate];

flux = cell(6,1);

% parameters
K1ac = p(1);    % or 0.02
K3fdp = p(2);
L3fdp = p(3)*1e6;
K3pep = p(4);
K2pep = p(5);
vemax = p(6);        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = p(7);        % or 0.45
ne = p(8);             % or 2
d = p(9);
V4max = p(10);
k1cat = p(11);   
V3max = p(12);    
V2max = p(13);   
ac = p(17);

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
ode = casadi.Function('ode',{x,p_all,ident_c,p_useless,acetate},{oderhs}); 

dfx_sym = [];
D2FX = [];





