% simulate ode using either SUNDIALS or ode45
function solution = run_ode(odefun, x0, fun_p)
if nargin<3
    fun_p = [];
end

if ~isempty(fun_p)
    odefunh = @(t,x)odefun(t,x,fun_p);
else
    odefunh = odefun;
end

%call to ode45 to simulate ode
[tout, yout] = ode45(odefunh, tspan, x0, options);

