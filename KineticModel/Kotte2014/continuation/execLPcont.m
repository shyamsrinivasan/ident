function [LPdata_f,LPdata_b] = execLPcont(xeq,apLP,ap,pvec,funame,data)
% perform orward and backward continuation from LPs
% data - data from equilibrium continuation
nvar = size(xeq,1);

% extract all limit points from data
id = cat(1,data.s1.index);
labels = cellstr(cat(1,data.s1.labels));
LPid = find(cellfun(@(x)strcmpi(x,'LP'),labels));

% get all limit points and corresponding parameters
xLP = data.x1(1:nvar,LPid);
pLP = data.x1(end,LPid);

% set continuation options
opt=contset;
opt=contset(opt,'MaxNumPoints',300);
opt=contset(opt,'MinStepSize',0.00001);
opt=contset(opt,'MaxStepSize',0.01);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Eigenvalues',1);

for ilp = 1:length(LPid)	
	pvec(apLP) = pLP(ilp);
	% change funame if you have to
	% newfname = @(x)funame(x,pvec);
	% continue forward from each limit point in data	
	LPdata_f = LPcontinuation(xLP(:,ilp),ap,funame,pvec,opt);
	% continue backward rom each point
	opt = contset(opt,'Backward',1);
	LPdata_b = LPcontinuation(xLP(:,ilp),ap,funame,pvec,opt);
end
