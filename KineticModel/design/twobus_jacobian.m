function jac = twobus_jacobian(x,p)
% x = (alpha,V)
% p = (P,Q);
jac = [-4*x(2)*cos(x(1)) -4*sin(x(1));...
        -4*x(2)*sin(x(1)) -8*x(2)+4*cos(x(1))];