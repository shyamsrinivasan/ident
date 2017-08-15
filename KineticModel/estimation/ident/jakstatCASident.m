function [ode,dfx_sym,D2FX,oderhs,x,p_all,ident_c,p_useless,acetate] =...
        jakstatCASident(idx)

x = casadi.SX.sym('x',4,1);
y = casadi.SX.sym('y',2,1);
tau = casadi.SX.sym('tau',1,1);
p_all = casadi.SX.sym('p_all',4,1);
ident_c = casadi.SX.sym('ident_c',1,1);
u = casadi.SX.sym('u',1,1);
s = casadi.SX.sym('s',2,1);

p = [p_all(1:idx-1);ident_c;p_all(idx:end)];

ode_state = [-p(1).*x(1).*u+2.*p(4).*x(4).^tau;...
            p(1).*x(1).*u-p(2).*x(2).^2;...
            p(2)./2.*x(2).^2-p(3).*x(3);...
            p(3).*x(3)-p(4).*x(4).^tau];

ode_output = [s(1).*(x(2)+2.*x(3));...
            s(2).*(x(1)+x(2)+2.*x(3))];
        
        