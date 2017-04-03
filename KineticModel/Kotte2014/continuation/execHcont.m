function Hdata = execHcont(funame,xeq,pvec,apH,ap,data,contpts)

if nargin<7
    contpts = 300;
end
% perform orward and backward continuation from LPs
% data - data from equilibrium continuation
nvar = size(xeq,1);

% extract all limit points from data
id = cat(1,data.s1.index);
label = cellstr(cat(1,data.s1.label));
Hid = id(cellfun(@(x)strcmpi(x,'H'),label));

% get all limit points and corresponding parameters
xH = data.x1(1:nvar,Hid);
pH = data.x1(end,Hid);

% set continuation options
opt=contset;
opt=contset(opt,'MaxNumPoints',contpts);
opt=contset(opt,'MinStepSize',0.00001);
opt=contset(opt,'MaxStepSize',0.1);
opt=contset(opt,'Singularities',1);
opt = contset(opt,'Eigenvalues',1);

% Hdata = cell(length(Hid),1);
% LPdata_b = cell(length(LPid),1);
ih = 1;
% for ilp = 1:length(LPid)	
    fprintf('\nHopf Point Continuation......\n');
	pvec(ap) = pH(ih);
	% change funame if you have to
	% newfname = @(x)funame(x,pvec);
	% continue forward from each limit point in data	
	Hdata = Hcontinuation(xH(:,ih),apH,funame,pvec,opt);
    Hdata.pH = pH(ih);
    Hdata.ap = ap;
    fprintf('Hopf Point Continuation Complete\n');
% end