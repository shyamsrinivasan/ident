function out = LorrenzMATCONT
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

function dM = fun_eval(t,kmrgd,sigma,rho,beta)
                     
pvec = [sigma,rho,beta];
dM = lorrenz(0,kmrgd,pvec);
    
% flux = KotteMATCONTflux(kmrgd,pvec);
% dM = zeros(3,1);
% % differential equations
% % PEP
% dM(1) = flux(1) - flux(4) - flux(5);
% % FBP
% dM(2) = flux(4) - flux(3);
% % enzymes
% % E
% dM(3) = flux(2) - d*kmrgd(3);


function [tspan,y0,options] = init
handles = feval(LorrenzMATCONT);

% obtain initial steady states
M = zeros(3,1);
M(1)  = 1;      % E
M(2)  = 1;   % PEP
M(3)  = 1;   % FBP

% substitute this with SUNDIALS
lorenzcall = @(t,x)lorrenz(t,x,[10;28;8/3]);
[~,yout] = ode45(lorenzcall,0:0.1:30,M);
y0 = yout(end,:);
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),...
                 'Hessians',handles(4),'HessiansP',handles(5));
tspan = [0 10];

function jac = jacobian(t,kmrgd,sigma,rho,beta)

function jacp = jacobianp(t,kmrgd,sigma,rho,beta)      
                     
function jacp = hessians(t,kmrgd,sigma,rho,beta)           
                     
function jacp = hessiansp(t,kmrgd,sigma,rho,beta)      

function jacp = der3(t,kmrgd,sigma,rho,beta)    
                     
function jacp = der4(t,kmrgd,sigma,rho,beta)    
                     
function jacp = der5(t,kmrgd,sigma,rho,beta)                         




