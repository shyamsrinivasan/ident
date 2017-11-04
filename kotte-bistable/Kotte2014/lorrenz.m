function dx = lorrenz(t,x,p)
dx = zeros(3,1);
sigma = p(1);
rho = p(2);
beta = p(3);

dx(1) = sigma*(x(2)-x(1));
dx(2) = rho*x(1)-x(2)-x(1)*x(3);
dx(3) = x(1)*x(2)-beta*x(3);