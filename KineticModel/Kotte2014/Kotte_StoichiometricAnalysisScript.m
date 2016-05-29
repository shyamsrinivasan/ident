%% analysis of SV = 0 for Kotte Model - 
% EM analysis using Cell Net Analyzer (CNA)

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

[hfig,hsubfig,prxnid,flag] = FluxEnvelope(FBAmodel,...
                        {'bmt2r','PEPt2r';...
                        'ACpts','ENZtr';...
                        'bmt2r','ACpts';...
                        'ACpts','PEPt2r'},...
                        {'ACt2r','ENZ1ex'});

[allefm,yield] = printEVyieldspace(hfig,hsubfig,{'bmt2r','PEPt2r';...
                        'ACpts','ENZtr';...
                        'EC_Biomass','ACpts';...
                        'ACpts','PEPt2r'},efm,idx,FBAmodel,find(strcmpi(model.rxns,'ACpts')));

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

% cut set calculations script
% getCNAMCS

