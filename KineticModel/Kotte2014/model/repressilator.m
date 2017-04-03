function out = repressilator
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

function dM = fun_eval(t,x,alpha,beta,delta,h)
% pi = [alpha,beta];
% ps = [delta,h];
p = [alpha,beta,delta,h];
dM = repressilatorNLAE(x,[],p);

function [tspan,y0,options] = init
handles = feval(repressilator);

% obtain initial steady states
y0 = zeros(1,6);
% substitute this with SUNDIALS
% odefun = @(t,x)designODE(t,x,pi,ps);
% [~,yout] = ode45(odefun,0:0.1:30,M);
% y0 = yout(end,:);
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),...
                 'Hessians',handles(4),'HessiansP',handles(5));
tspan = [0 10];

function jac = jacobian(t,x,alpha,beta,delta,h)
p = [alpha,beta,delta,h];
jac = getKotteJacobian(@repressilatorNLAE,x,p,[]);

function jacp = jacobianp(t,x,alpha,beta,delta,h)      
                     
function jacp = hessians(t,x,alpha,beta,delta,h)           
                     
function jacp = hessiansp(t,x,alpha,beta,delta,h)      

function jacp = der3(t,x,alpha,beta,delta,h)    
                     
function jacp = der4(t,x,alpha,beta,delta,h)    
                     
function jacp = der5(t,x,alpha,beta,delta,h)          