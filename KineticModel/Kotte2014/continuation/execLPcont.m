function [LPdata_f,LPdata_b] = execLPcont(funame,xeq,pvec,apLP,ap,data,contpts)
if nargin<7
    contpts = 300;
end
% perform orward and backward continuation from LPs
% data - data from equilibrium continuation
nvar = size(xeq,1);

% extract all limit points from data
id = cat(1,data.s1.index);
label = cellstr(cat(1,data.s1.label));
LPid = id(cellfun(@(x)strcmpi(x,'LP'),label));

% get all limit points and corresponding parameters
xLP = data.x1(1:nvar,LPid);
pLP = data.x1(end,LPid);

% set continuation options
opt=contset;
opt=contset(opt,'MaxNumPoints',contpts);
opt=contset(opt,'MinStepSize',0.00001);
opt=contset(opt,'MaxStepSize',0.01);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Eigenvalues',1);

LPdata_f = cell(length(LPid),1);
LPdata_b = cell(length(LPid),1);

for ilp = 1:length(LPid)	
    fprintf('\nLimit Point Continuation #%d of %d......\n',ilp,length(LPid));
	pvec(ap) = pLP(ilp);
	% change funame if you have to
	% newfname = @(x)funame(x,pvec);
	% continue forward from each limit point in data	
	LPdata_f{ilp} = LPcontinuation(xLP(:,ilp),apLP,funame,pvec,opt);
    LPdata_f{ilp}.pLP = pLP(ilp);
    LPdata_f{ilp}.ap = ap;
	% continue backward rom each point
	opt = contset(opt,'Backward',1);
	LPdata_b{ilp} = LPcontinuation(xLP(:,ilp),apLP,funame,pvec,opt);
    LPdata_b{ilp}.pLP = pLP(ilp);
    LPdata_b{ilp}.ap = ap;
    fprintf('Limit Point Continuation #%d of %d Complete\n',ilp,length(LPid));
end
