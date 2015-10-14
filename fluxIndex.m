% function [Vind,Vuptake,VFup,VFex,Vex,bmrxn,Vup,Vdn] = fluxIndex(model,nt_rxn,newS)
function [Vind,VFex,Vex,bmrxn] = fluxIndex(model,nt_rxn,newS)
if nargin<3
    newS = model.S;
end
if nargin<2
    nt_rxn = size(newS,2);
end
%external mets
% exind = ~cellfun('isempty',regexp(model.mets,'\w(?:\[e\])$'));

[~,allrxns] = find(newS(:,1:nt_rxn));
all_unqrxns = unique(allrxns);
nmetab_allrxns = histc(allrxns,all_unqrxns);
ex_rxn = all_unqrxns(nmetab_allrxns == 1);%one sided rxns

%FBA style Reactions A[e/c] <==> | <==> A[e/c]
%All one sided reactions
[~,rxn] = find(newS(:,ex_rxn)~=0);
VFex = ex_rxn(rxn);

% %one sided rxns
% [~,rxn] = find(newS(:,ex_rxn)<0);
% VFex = ex_rxn(rxn);%Excrretion from the system
% [~,rxn] = find(newS(:,ex_rxn)>0);
% VFup = ex_rxn(rxn);%Uptake into the system
% try
%     VFext = [VFup VFex];
% catch
%     VFext = [VFup;VFex];
% end

%All exchange reactions: A[e] ---> A[c] | A[e] ---> B[c] |
%A[e] + B[c] ---> A[c] + D[c] 
%and Matching VFex with Vex for met[e]
Vex = [];
for ivf = 1:length(VFex)
    %Find met[e]
    metid = logical(newS(:,VFex(ivf))~=0);
    %Find Vex for met[e]
    noex_rxn = setdiff(1:nt_rxn,VFex(ivf));
    [~,rxn] = find(newS(metid,noex_rxn)~=0);
    Vex = union(Vex,noex_rxn(rxn));
end

% noex_rxn = setdiff(1:nt_rxn,ex_rxn);
% [~,rxn] = find(newS(exind,noex_rxn)~=0);
% Vex = noex_rxn(rxn);

% %Uptake Reactions:
% noex_rxn = setdiff(1:nt_rxn,ex_rxn);
% [~,rxn] = find(newS(exind,noex_rxn)<0);
% Vuptake = noex_rxn(rxn);
% 
% %Excrertion Reaction: A[c] ---> A[e] | A[c] ---> B[e] | A[c] + B[c] --->
% %A[e] + D[c]
% [~,rxn] = find(newS(exind,noex_rxn)>0);
% Vex = noex_rxn(rxn);




%Matching VFup with Vuptake
% Vup = [];
% for ivf = 1:length(VFup)
%     nmet = newS(:,VFup(ivf))>0;
%     [~,rxn] = find(newS(nmet,noex_rxn)<0);
%     Vup = union(Vup,noex_rxn(rxn));
% end


%Remove Vex from considering A[c] + B[e] ---> C[c] + D[c] type reactions
% remVex = zeros(length(Vex),1);
% for iex = 1:length(Vex)
%     if length(find(newS(:,Vex(iex))>0))>1 ||... 
%        length(find(newS(:,Vex(iex))<0))>1
%         remVex(iex) = 1;        
%     end
% end
% Vex(logical(remVex)) = [];


%matchinf VFex with Vex
% Vdn = [];
% for ive = 1:length(VFex)
%     nmet = newS(:,VFex(ive))<0;
%     [~,rxn] = find(newS(nmet,noex_rxn)>0);
%     if length(rxn)>1
%         remVdn = zeros(length(rxn),1);
%         for ir = 1:length(rxn)
%             if length(find(newS(:,noex_rxn(rxn(ir)))>0))>1 ||...
%                length(find(newS(:,noex_rxn(rxn(ir)))<0))>1
%                remVdn(ir) = 1;
%             end
%         end
%         rxn(logical(remVdn)) = [];
%     end
%     Vdn = d_uUnion(Vdn,noex_rxn(rxn));
% end

%Identify biomass reaction
% [~,bmrxn] = find(newS(strcmpi(model.mets,'Biomass'),:) > 0);
bmrxn = find(strcmpi(model.rxns,'EC_biomass'));
VFex = setdiff(VFex,bmrxn);
Vex = setdiff(Vex,bmrxn);

try
%     Vind = setdiff(1:nt_rxn,[Vuptake;bmrxn;VFext;Vex']);%intracellular rxns
    Vind = setdiff(1:nt_rxn,[VFex;bmrxn;Vex']);%intracellular rxns
catch
%     Vind = setdiff(1:nt_rxn,[Vuptake bmrxn VFext' Vex]);%intracellular rxns
    Vind = setdiff(1:nt_rxn,[VFex' bmrxn Vex']);%intracellular rxns
    
end
end