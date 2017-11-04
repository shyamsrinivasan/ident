function fx = bifurequations(x)
xstate = x(1:nvar);
x_w1 = x(nvar+1:2*nvar);
x_l1 = x(2*nvar+1:end);

newlambda = oldlambda + x_l1*n;
fx_state = twobusAE(xstate,newlambda);
fx_jacobian = x_w1'*jacobian(xtate,newlambda);
fx_nzw1 = w1'*c - 1;

fx = [fx_state;fx_jacobian;fx_nzw1];