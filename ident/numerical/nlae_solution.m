function solution = nlae_solution(fun, p, initial_value, typical_vals)
if nargin<4
    typical_vals = [];
end
if nargin<3
    initial_value = ones(3, 1);
end
if ~isempty(p)
    nlaefun = @(x)fun(x, p);
else
    nlaefun = fun;
end

% set fsolve options
options = optimoptions('fsolve','MaxIter',500,...
                                'TolFun',1e-12,...
                                'TolX',1e-12);
if ~isempty(typical_vals)
    options = optimoptions(options,'TypicalX',typical_vals,...
                                    'Display','iter',...
                                    'MaxIter',9000,...
                                    'MaxFunEvals',20000);
end
% use fsolve to solve NLAE equation for parameters for flux v3
[x, fval, exitflag, output] = fsolve(nlaefun, initial_value, options);

if exitflag
    solution.p = x;
    solution.res = fval;
    solution.message = output.message;
else
    solution.p = [];
    solution.res = [];
    solution.message = [];
end
solution.flag = exitflag;

return
