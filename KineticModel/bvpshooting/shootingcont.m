% continuation for 2 point BVP
% define problem
% Holt.m
% initial conditions
nvar = 5;
yinit = zeros(nvar,1);

% known initial conditios @ t0
yiknwn = zeros(nvar,1);
yiknwn(1) = 1; yiknwn(2) = 1; yiknwn(4) = 1;
% known terminal conditions @ tf
yterm = zeros(nvar,1);
yterm(4) = 1;

yiknwn = find(yiknwn);
r = length(yiknwn);
yiunkwn = setdiff(1:nvar,yiknwn);

% terminal conditions @ tf
yfknwn = zeros(nvar,1);
yfknwn(2) = 1; yfknwn(4) = 1;
yfunkwn = find(~yfknwn);
yfknwn = find(yfknwn);

% choose unknown initial conditions
yinit(yiunkwn) = [-.978197694;0.646786696];
delyi = getvaldiff(yinit,yinit);

% cotinuation parameters
F = 1.0; 
tau = 0.2;
ti = 0;
tf = 3.5;
tterm = 11.0;
eps = 1e-4;

% integrate till time t1 till which there are no numerical problems
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
% [yi,yf,delyi,delyf] =...
% itershooting(@HoltODE,yinit,yterm,ti,tf,yiunkwn,yfknwn,delyi,[],[],opts);
yi = yinit;
[tf,yf] = integrateshooting(@HoltODE,ti,tf,yi,opts);
% [~,ydyn] = ode45(@HoltODE,ti:0.1:tf,yi,opts);
% yf = ydyn(end,:)';
delyf = getvaldiff(yterm,yf);

allyi = yi;
allti = ti;
alltf = tf;
while tf<tterm
    flag = 1;
    while flag        
        % solve 2 point bvp over (t0,t1) using goodman lance method
        [yi,yf,tf,delyi,~,flag] =...
        execshooting(@HoltODE,yi,yterm,ti,tf,yiunkwn,yfknwn,delyi,delyf,yf,eps);    
        if flag
            % store yi,ti,tf values
            allyi = [allyi yi];
            allti = [allti;ti];
            alltf = [alltf;tf];
            tf = tf + tau;
        end     
        yfold = yf;
        [~,yf] = integrateshooting(@HoltODE,ti,tf,yi,opts);        
        delyf = getvaldiff(yterm,yf);
%         [yi,yf,delyi,delyf] =...
%         itershooting(@HoltODE,yi,yf,ti,tf,yiunkwn,yfknwn,delyi,newdelyf,ysimf,opts);
    end
    tau = tau/2;
end





% extend time interval from (t0,t1) to (t0,t2)