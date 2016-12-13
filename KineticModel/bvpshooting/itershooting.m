function [yinew,ysimf,delyf] =...
itershooting(fh,yi,yf,ti,tf,yiunkwn,yfknwn,delyi,delyf,ysimf,opts,r,nvar)
% perform kth iteration of the goodman-lance shooting method to get delyf -
% the difference in terminal boundary conditions
if nargin<13
    nvar = size(yi,1);
end
if nargin<12
    yiknwn = setdiff(1:nvar,yiunkwn);
    r = length(yiknwn);
end
if nargin<11
    opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
end
if nargin<10 || isempty(ysimf)
    [~,ydyn] = ode45(fh,ti:0.1:tf,yi,opts);
    ysimf = ydyn(end,:)';
end    
if nargin<9 ||isempty(delyf)
    delyf = zeros(size(yi,1),1);
end


% integrate adjoint equation in reverse time from tf to t0 for m =
% 1,...,n-r
% fprintf('Solving adjoint Equation...\n'); 
xic = zeros(nvar,nvar-r);
for m = 1:nvar-r
    getAdj = @(t,x)adjointEquation(t,x,ysimf);
    xf = zeros(nvar,1);
    xf(yfknwn(m)) = 1;
    [~,xint] = ode45(getAdj,ti:-0.1:-tf,xf,opts);
    xic(:,m) = xint(end,:)';
end
    
% solve n-r algebraic equations given by 
% fprintf('Solving algebraic Equations...\n');
delyvar = delyi(yiunkwn);
getAgeq = @(delyi)BVPshooting_algebraic(delyi,delyf,xic,yiunkwn,yfknwn);
fx = getAgeq(delyvar);
options = optimoptions('fsolve','Display','off',...
                       'TolFun',1e-10,...
                       'TolX',1e-10,...
                       'MaxFunEvals',1000000,...
                       'MaxIter',50000);
[new_delyi,fval,exitflag,output,jacobian] = fsolve(getAgeq,delyvar,options);
delyi(yiunkwn) = new_delyi(:);
yinew = yi + delyi;
    
% integrate system with guessed/new intial conditions
% fprintf('Integrated system to find final value...\n');
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
[~,ydyn] = ode45(fh,ti:0.1:tf,yinew,opts);
ysimf = ydyn(end,:)';
    
% check difference in final values
delyf = getvaldiff(yf,ysimf); % assuming ysimf is obtained at tf