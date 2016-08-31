% getCNAMCS
% get MCS, cMCS and regCMCS - constrained minimal cut sets - 
% ends with an  error due to lack of MCSs
cnap.reacMin(cnap.reacMin == -100) = -1; % Inf;
cnap.reacMax(cnap.reacMax == 100) = 1; % Inf;
cnap.reacMax(strcmpi(cellreacID,'ACt2r')) = 0;

tar = zeros(1,cnap.numr);
nT = 1;
T = repmat(tar,nT,1);
% 
cellreacID = cellstr(cnap.reacID);
cellreacID = cellfun(@(x)strtrim(x),cellreacID,'UniformOutput',false);
T(1,[find(strcmpi(cellreacID,'PEPt2r')) find(strcmpi(cellreacID,'ACt2r'))]) =...
  [1 0.1];
% T(1,strcmpi(cellreacID,'EC_Biomass')) = -1;
% T(3,strcmpi(cellreacID,'ACpts')) = -1;
t = zeros(nT,1);
% t(1) = 0.5;
% t(1) = -0.5;
% t(3) = -1;
t(1) = 0;

nD = 1;
D = repmat(tar,nD,1);
D(1,strcmpi(cellreacID,'bmt2r')) = -1;
% D(2,strcmpi(cellreacID,'PEPt2r')) = -1;
% D(2,strcmpi(cellreacID,'ACt2r')) = 1;
% D(1,strcmpi(cellreacID,'bmt2r')) = -1;
d = zeros(nD,1);
d(1) = -0.7;
% d(2) = -0.19;
% d(3) = 0;
% d(1) = -0.5;

notknockable = [1 2 3 5 6 9];
maxMCS = 100;
maxMCSsize = 4;
filename = [];

% mcs =...
% CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCS,maxMCSsize,filename);
% printCS(mcs,cnap.reacID);

%regMCS
regulation.reg_ind = [1 2 3];
regulation.reg_down_up = ones(2,length(regulation.reg_ind));
regulation.numregsteps = 3;
[mcs,reacNames] = CNAregMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCS,maxMCSsize,...
                          filename,0,[],regulation);
