function [lb,ub] = ident_bounds(np,lb,ub)
if nargin<3
    ub = zeros(np,1);
end
if nargin<2
    lb = zeros(np,1);
end

% k1cat - PLE
lb(:) = .05;
lb(3) = 2;
lb(6) = .05;
lb(7) = .1;
lb(8) = .1;

ub(:) = 2;
ub(3) = 5;
ub(6) = 1;
