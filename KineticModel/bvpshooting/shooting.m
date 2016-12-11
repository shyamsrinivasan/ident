% continuation with shooting methods for 2 point BVPs
% 1. Agarwal, R., 1979
% 2. Shipman, J. S., and Roberts, S. M., 1967

% shooting method of Lance and Goodman
% variables
% nvar     - number of variables
% yiknwn    - rx1 variables for which initial conditions are given
% yfknwn    - (nvar-r)x1 variables for which terminal conditions are given
% yinit     - initial values (guesses and given)
% yterm     - terminal values (given values only)
% dely      - [delyi delyf]
% delyi     - difference between given and calculate initial conditions
%           - is 0 initially when strating at the first iteration
% delyf     - difference between given and calculated terminal conditions
% xic       - nx1 adjoint variables obtained from reverse integration of
%             adjoint equation
% xf        - nx1 adjoint variable initial conditions for reverse
%             integration


% define problem
% Holt.m
% initial conditions
nvar = 5;
yinit = zeros(nvar,1);

% known initial conditios @ t0
yiknwn = zeros(nvar,1);
yiknwn(1) = 1; yiknwn(2) = 1; yiknwn(4) = 1;
yiknwn = find(yiknwn);
r = length(yiknwn);
yiunkwn = setdiff(1:nvar,yiknwn);

% terminal conditions
yterm = zeros(nvar,1);
yterm(4) = 1;

% known terminal conditions @ tf
yfknwn = zeros(nvar,1);
yfknwn(2) = 1; yfknwn(4) = 1;
yfunkwn = find(~yfknwn);
yfknwn = find(yfknwn);

% guess unknown initial conditions
yinit(yiunkwn) = [-1;0.6];

% beginning of iterative loop
% assume dely(t0) = 0
delyi = getvaldiff(yinit,yinit);

% integrate system with guessed intial conditions
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
[tf,ydyn] = ode45(@HoltODE,0:0.1:3.5,yinit,opts);
yf = ydyn(end,:)';

% check difference in final values
delyf = getvaldiff(yterm,yf); % assuming yf is obtained at tf

% integrate adjoint equation in reverse time from tf to t0 for m =
% 1,...,n-r
xic = zeros(nvar,nvar-r);
for m = 1:nvar-r
    getAdj = @(t,x)adjointEquation(t,x,yf);
    xf = zeros(nvar,1);
    xf(yfknwn(m)) = 1;
    [~,xint] = ode45(getAdj,0:-0.1:-3.5,xf,opts);
    xic(:,m) = xint(end,:)';
end

% solve n-r algebraic equations given by 
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
yinew = yinit + delyi;






