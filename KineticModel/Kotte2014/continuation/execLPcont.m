function execLPcont(xeq,ap,pvec,funame,data)
% perform orward and backward continuation from LPs
% data - data from equilibrium continuation

% extract all limit points from data
% continue from each limit point in data
data = LPcontinuation(xeq,ap,funame,pvec)
