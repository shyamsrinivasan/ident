function [yinew,delyi] =...
itershooting(fh,ti,tf,yi,yf,delyi,delyf,yiunkwn,yfknwn,nvar,r,opts)
% perform kth iteration of the goodman-lance shooting method to get delyf -
% the difference in terminal boundary conditions
if nargin<12
    opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
end
% if nargin<13
%     nvar = size(yi,1);
% end
% if nargin<12
%     yiknwn = setdiff(1:nvar,yiunkwn);
%     r = length(yiknwn);
% end
% if nargin<11
%     opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
% end
% if nargin<10 || isempty(ysimf)
%     [~,ydyn] = ode45(fh,ti:0.1:tf,yi,opts);
%     ysimf = ydyn(end,:)';
% end    
% if nargin<9 ||isempty(delyf)
% %     delyf = zeros(size(yi,1),1);
%     delyf = getvaldiff(yf,ysimf);
% end

% reverse integrate adjoint equations
xic = zeros(nvar,nvar-r);
for m = 1:nvar-r
    getAdj = @(t,x)adjointEquation(t,x,yf);
    xf = zeros(nvar,1);
    xf(yfknwn(m)) = 1;
    [tint,xint] = ode15s(getAdj,tf:-0.1:ti,xf,opts);
    xic(:,m) = vpa(xint(end,:)',16);
end

% solve n-r algebraic equations
delyvar = delyi(yiunkwn);
getAgeq = @(delyi)BVPshooting_algebraic(delyi,delyf,xic,yiunkwn,yfknwn);
fx = getAgeq(delyvar);
options = optimoptions('fsolve','Display','off',...
                       'TolFun',1e-16,...
                       'TolX',1e-16,...
                       'MaxFunEvals',1000000,...
                       'MaxIter',50000);
new_delyi = fsolve(getAgeq,delyvar,options);

delyi(yiunkwn) = new_delyi(:);

% get new yi
yinew = yi + delyi;


% % integrate adjoint equation in reverse time from tf to t0 for m =
% % 1,...,n-r
% % fprintf('Solving adjoint Equation...\n'); 
% xic = zeros(nvar,nvar-r);
% for m = 1:nvar-r
%     getAdj = @(t,x)adjointEquation(t,x,ysimf);
%     xf = zeros(nvar,1);
%     xf(yfknwn(m)) = 1;
%     [tint,xint] = ode15s(getAdj,tf:-0.1:ti,xf,opts);
%     xic(:,m) = vpa(xint(end,:)',16);
% end
%     
% % solve n-r algebraic equations given by 
% % fprintf('Solving algebraic Equations...\n');
% delyvar = delyi(yiunkwn);
% getAgeq = @(delyi)BVPshooting_algebraic(delyi,delyf,xic,yiunkwn,yfknwn);
% fx = getAgeq(delyvar);
% options = optimoptions('fsolve','Display','off',...
%                        'TolFun',1e-16,...
%                        'TolX',1e-16,...
%                        'MaxFunEvals',1000000,...
%                        'MaxIter',50000);
% [new_delyi,fval,exitflag,output,jacobian] = fsolve(getAgeq,delyvar,options);
% delyi(yiunkwn) = vpa(new_delyi(:),16);
% yinew = yi + delyi;
    
% % integrate system with guessed/new intial conditions
% % fprintf('Integrated system to find final value...\n');
% [tf,ysimf,status] = integrateshooting(ti,tf,yi,opts);
    
% % check difference in final values
% delyf = getvaldiff(yf,ysimf); % assuming ysimf is obtained at tf