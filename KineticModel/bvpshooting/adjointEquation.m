function dx = adjointEquation(t,x,y)

% get jacobian
J = HoltJacobian(y);

% get adjoint equation
% dx(1) = -[dg1/dy1*x1 + dg2/dy1*x2 + ... + dgn/dy1*xn]
dx = -(x*J);
dx = dx';