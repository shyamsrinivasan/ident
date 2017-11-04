function [cnap,errval] = getCNAmodel(FBAmodel)
% get CNA model
% run after running Kotte2014_script or from it
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