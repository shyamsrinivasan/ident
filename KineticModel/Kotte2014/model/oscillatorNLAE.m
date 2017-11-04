function dX = oscillatorNLAE(x,model,p)
dX = zeros(3,1);
dX = cons(dX,x); % for use with ADMAT for jacobian calculation

q1 = p(1);
q2 = p(2);
q3 = p(3);
q4 = p(4);
q5 = p(5);
q6 = p(6);
k = p(7);

z = 1-x(1)-x(2)-x(3);
dX(1) = 2*q1*z^2-2*q5*x(1)^2-q3*x(1)*x(2);
dX(2) = q2*z-q6*x(2)-q3*x(1)*x(2);
dX(3) = q4*z-k*q4*x(3);