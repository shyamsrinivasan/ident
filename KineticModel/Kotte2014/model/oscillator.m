function out = oscillator
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

function dM = fun_eval(t,kmrgd,q1,q2,q3,q4,q5,q6,k)
p = [q1;q2;q3;q4;q5;q6;k];
dM = oscillatorNLAE(kmrgd,p);

function [tspan,y0,options] = init
handles = feval(oscillator);

% obtain initial steady states
y0 = zeros(1,3);
% substitute this with SUNDIALS
% odefun = @(t,x)oscillatorODE(t,x,p);
% [~,yout] = ode45(odefun,0:0.1:100,M);
% y0 = yout(end,:)';
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),...
                 'Hessians',handles(4),'HessiansP',handles(5));
tspan = [0 10];

function jac = jacobian(t,x,q1,q2,q3,q4,q5,q6,k)

function jacp = jacobianp(t,x,q1,q2,q3,q4,q5,q6,k)      
                     
function jacp = hessians(t,x,q1,q2,q3,q4,q5,q6,k)           
                     
function jacp = hessiansp(t,x,q1,q2,q3,q4,q5,q6,k)      

function jacp = der3(t,x,q1,q2,q3,q4,q5,q6,k)    
                     
function jacp = der4(t,x,q1,q2,q3,q4,q5,q6,k)    
                     
function jacp = der5(t,x,q1,q2,q3,q4,q5,q6,k)          