function [p_bounds_lb,p_bounds_ub] = p_bounds()

% all parameters available (ordered)
plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep','vemax','KeFDP','ne',...
        'd','V4max','k1cat','V3max','V2max'}; 
    
p_bounds_lb = zeros(length(plist),1);
p_bounds_ub = zeros(length(plist),1);
% Kms
p_bounds_lb([1 2 4 5 7]) = .008;
p_bounds_ub([1 2 4 5 7]) = 1;
% kcats/Vmax
p_bounds_lb([6 10 11 12 13]) = .1;
p_bounds_ub([6 10 11 12 13]) = 2;
% misc
p_bounds_lb([3 8 9]) = [1;1;0];
p_bounds_ub([3 8 9]) = [5;3;1];
