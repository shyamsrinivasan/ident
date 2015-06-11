% function [newmodel] = reduceModelSize(oldmodel,regprotein,reglist)
% Reduce model to include only specific regulators/genes
function [newmodel] = reduceModelSize(oldmodel,regprotein,reglist,genelist)
if nargin < 4
    genelist = {};
end

if nargin < 3%??
    fprintf('No list to parse through\n');
    newmodel = oldmodel;
    return
end

%Remove all genes other than those in the reglist/Remove all genes on the
%list??
%======================================================================
if ~isempty(genelist)
    
end     
%     %Determine genes regulated by TFs ArcA and FNR and remove them
%     %====================================================================
%     rtf = strcmp(reglist{ireg},trnmodel.Protein);
%     %gindx = logical(trnmodel.RS(:,rtf));
%     newmodel.RS(logical(trnmodel.RS(:,rtf)),rtf) = 0;

%Remove all proteins other than those in the reglist
%======================================================================
if ~isempty(reglist)
    % remove_prot = setdiff(trnmodel.Protein,reglist);
    remove_prot = setdiff(regprotein,reglist);
    nremove_prot = length(remove_prot);

    newmodel = blacklist_tfs(oldmodel,remove_prot,regprotein);

    prot_col = zeros(length(oldmodel.Protein),1);
    gaprot_col = zeros(length(regprotein),1);

    for iremove = 1:nremove_prot
        rpindx = strcmp(remove_prot{iremove},newmodel.Protein);
        rpindx1 = strcmp(remove_prot{iremove},regprotein);
        newmodel.Protein{rpindx} = {};
        prot_col(rpindx) = 1;
        gaprot_col(rpindx1) = 1;
    end
   
    prot_col = logical(prot_col);
    gaprot_col = logical(gaprot_col);

    newmodel.Protein = newmodel.Protein(~cellfun('isempty',newmodel.Protein));
    %newmodel.Protein(logical(prot_col)) = {};
    newmodel.RS(:,prot_col) = [];
    newmodel.Coefficient(:,prot_col) = [];
    newmodel.Kmax(prot_col,:) = [];
    newmodel.Ks(prot_col,:) = [];

    newmodel.trate(gaprot_col,:) = [];
end

%fprintf('End of function\n');
end