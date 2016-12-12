function delyf = itershooting(fh,yi,yf,yiunkwn,yfknwn,delyi,delyf,ysimf,opts)
if nargin<9
    opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
end

if nargin<8
    [~,ydyn] = ode45(fh,0:0.1:3.5,yi,opts);
    ysimf = ydyn(end,:)';
end
    

% integrate adjoint equation in reverse time from tf to t0 for m =
% 1,...,n-r
fprintf('Solving adjoint Equation...\n'); 
xic = zeros(nvar,nvar-r);
for m = 1:nvar-r
    getAdj = @(t,x)adjointEquation(t,x,ysimf);
    xf = zeros(nvar,1);
    xf(yfknwn(m)) = 1;
    [~,xint] = ode45(getAdj,0:-0.1:-3.5,xf,opts);
    xic(:,m) = xint(end,:)';
end
    
% solve n-r algebraic equations given by 
fprintf('Solving algebraic Equations...\n');
delyvar = delyi(yiunkwn);
getAgeq = @(delyi)BVPshooting_algebraic(delyi,delyf,xic,yiunkwn,yfknwn);
fx = getAgeq(delyvar);
options = optimoptions('fsolve','Display','iter',...
                       'TolFun',1e-10,...
                       'TolX',1e-10,...
                       'MaxFunEvals',1000000,...
                       'MaxIter',50000);
[new_delyi,fval,exitflag,output,jacobian] = fsolve(getAgeq,delyvar,options);
delyi(yiunkwn) = new_delyi(:);
yi = yi + delyi;
    
% integrate system with guessed/new intial conditions
fprintf('Integrated system to find final value...\n');
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
[~,ydyn] = ode45(fh,0:0.1:3.5,yi,opts);
ysimf = ydyn(end,:)';
    
% check difference in final values
delyf = getvaldiff(yf,ysimf); % assuming ysimf is obtained at tf