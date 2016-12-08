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
r = length(yknwn);

% terminal conditions
yterm = zeros(nvar,1);
yterm(4) = 1;

% known terminal conditions @ tf
yfknwn = zeros(nvar,1);
yfknwn(2) = 1; yfknwn(4) = 1;
yfknwn = find(yfknwn);

% guess unknown initial conditions
yinit(~uiknwn) = [-1;0.6];

% integrate system with guessed intial conditions
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
[tf,yf] = ode45(@Holt,0:0.1:3.5,yinit,opts);

% check difference in final values
delyf = yterm - yf; % assuming yf is obtained at tf

% integrate adjoint equation in reverse time from tf to t0 for m =
% 1,...,n-r
xic = zeros(nvar,nvar-r);
for m = 1:nvar-r
    getAdj = @(t,x)adjointEquation(t,x,y);
    xf = zeros(nvar,1);
    xf(yfknwn(m)) = 1;
    [~,xint] = ode45(getAdj,-3.5:-0.1:0,xf,opts);
    xic(:,m) = xint(:,end);
end

% solve n-r algebraic equations
delyf = zeros(nvar-r,1);
for m = 1:nvar-r
    delyf(m) = sum(xic(yfknwn,m).*delyi(yfknwn))
end






