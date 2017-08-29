function [lb,ub] = ident_bounds(np)

lb = zeros(np,1);
ub = zeros(np,1);

% k1cat - PLE
lb(:) = .01;
lb(3) = 1;
lb(10) = .05;
lb(11) = .1;
lb(12) = .1;

ub(1) = 2;
ub(2) = 2;
ub(3) = 5;
ub(4) = 2;
ub(5) = 2;
ub(6) = 2;
ub(7) = 4;
ub(8) = 4;
ub(9) = 1;
ub(10) = 1;
ub(11) = 2;
ub(12) = 2;
