%% analysis of SV = 0 for Kotte Model - 
% EM analysis using Cell Net Analyzer (CNA)

% convert model to conform to CNA form
spec = ones(size(FBAmodel.S,1),1)';
spec(1:FBAmodel.nint_metab) = 0;
cnap.has_gui = 0;
cnap.net_var_name = 'KotteGmodel';
cnap.type = 1;
cnap.specID = char(FBAmodel.mets); 
cnap.specLongName = char(FBAmodel.mets);
cnap.specExternal = spec;
cnap.specInternal = find(~cnap.specExternal);
cnap.nums = size(FBAmodel.S,1);
cnap.numis = size(cnap.specInternal,2);
% cnap.macroID = 'BC1';
% cnap.macroLongName = 'BC1';
% cnap.macroComposition =...
% sparse(find(strcmpi(FBAmodel.mets,'bm[c]')),1,1,length(FBAmodel.mets),1);
% cnap.macroDefault = 1;
% cnap.nummac = 1;
cnap.stoichMat = full(FBAmodel.S); % zeros(length(FBAmodel.mets),1)];
cnap.numr = size(cnap.stoichMat,2);
cnap.reacID = char(FBAmodel.rxns); % char([FBAmodel.rxns;'mue']);
cnap.objFunc = zeros(cnap.numr,1);
cnap.reacMin = FBAmodel.vl; % 0];
cnap.reacMax = FBAmodel.vu; % 100];

[cnap,errval] = CNAgenerateMFNetwork(cnap);

% save CNA model
% cnap.path =...
% 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\KotteGmodel';
% cnap = CNAsaveNetwork(cnap);

% flux optimization
% constr = zeros(cnap.numr,1);
% constr(constr==0) = NaN;
% constr(1) = 10;
% constr(4) = 10;
% constr(9) = 0;
% constr(14) = 1;
% [flux,success,status] = CNAoptimizeFlux(cnap,constr,cnap.macroDefault,2,2);

% flux variablity
% reacval = zeros(cnap.numr,1);
% reacval(reacval==0) = NaN;
% reacval(5) = -10;
% reacval(11) = -10;
% [minFlux,maxFlux,success,status] =...
% CNAfluxVariability(cnap,reacval,cnap.macroDefault,2);

% phase plane analysis
% status =
% CNAplotPhasePlane(cnap,reacval,cnap.macroDefault,[14;10;12;13],2);

%% EFM calculation
constr = zeros(cnap.numr,4);
constr(constr==0) = NaN;
% exclude efms w/ exchanges only
constr([5 10 11 12 13],1) = 0;
% include efms with only all major reactions
% constr([8],1) = 1;
% constr(1,2) = 1;
% constr(9,2) = 0.1961;
% constr(6,2) = 0.3;
constr(8,2) = 1;
constr(1,3) = 1;
% constr(7,2) = -1;
% % constr(14,2) = 0.1;
% constr(9,2) = 0.1;
% constr([5 7],3) = 0;
% constr([1 4],4) = 1;

[efm,rev,idx,ray] = CNAcomputeEFM(cnap,constr,2,1,0,0);
printEFM(efm,idx,ray,cnap);

[allefm,yield] = printEVyieldspace([],hsubfig,prxnid,efm,idx,FBAmodel,find(strcmpi(model.rxns,'ACpts')));

% cut set calculation
% reacID = cnap.reacID(idx,:);
% target = efm([1 4],:);
% set2save(1).tabl2save = efm([1 2],:);
% set2save(2).tabl2save = efm(2,:);
% set2save(1).min2save = 2;
% set2save(2).min2save = 1;
% set2save = [];
% cutsets = CNAcomputeCutsets(target,10,reacID,set2save);
% printCS(cutsets,reacID);

%% constrained minimal cut sets - 
% ends with an  error due to lack of MCSs
cnap.reacMin(cnap.reacMin == -100) = -Inf;
cnap.reacMax(cnap.reacMax == 100) = Inf;

tar = zeros(1,cnap.numr);
nT = 1;
T = repmat(tar,nT,1);
% 
cellreacID = cellstr(cnap.reacID);
cellreacID = cellfun(@(x)strtrim(x),cellreacID,'UniformOutput',false);
T(1,strcmpi(cellreacID,'PEPt2r')) = 1;
% T(3,strcmpi(cellreacID,'ACpts')) = -1;
t = zeros(nT,1);
% t(1) = 0.5;
t(1) = 0.19;
% t(3) = -1;

nD = 3;
D = repmat(tar,nD,1);
D(1,strcmpi(cellreacID,'ACpts')) = -1;
D(2,strcmpi(cellreacID,'PEPt2r')) = -1;
D(3,strcmpi(cellreacID,'bmt2r')) = 1;
d = zeros(nD,1);
d(1) = 0;
d(2) = -0.19;
d(3) = -0.2;

notknockable = [];
maxMCS = 100;
maxMCSsize = 5;
filename = [];

mcs =...
CNAMCSEnumerator(cnap,T,t,[],[],notknockable,maxMCS,maxMCSsize,filename);
printCS(mcs,cnap.reacID);
