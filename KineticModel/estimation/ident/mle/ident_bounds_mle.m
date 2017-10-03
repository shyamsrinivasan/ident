function [lb,ub] = ident_bounds_mle(np,lb,ub)
if nargin<3
    ub = zeros(np,1);
end
if nargin<2
    lb = zeros(np,1);
end

% k1cat - MLE holding transcription parameters and acetate constant 
% {'vemax','KeFDP','ne','d','acetate'}
lb(:) = .05;
lb(3) = 1;
lb(6) = .1;
lb(7) = .1;
lb(8) = .1;
lb(9) = .1;

ub(:) = 3;
ub(1) = 2;
ub(2) = 2;
ub(3) = 4.5;
ub(4) = 2;
ub(5) = 2;
ub(6) = 1;
ub(7:9) = 1.5;

% ub(10) = 2;
% ub(11) = 2;
% ub(13) = 1;
