%function [newmodel,status] = blacklist_tfs(model,list,blacklist,recurmodel)
%**************************************************************************
% Usage 1: Reduce model appicability to selected environmental conditions        
% [model] = blacklist_tfs(fullmodel,metab,blacklist_flag) 

% Input
% fullmodel - Model to be reduced from
% metab - List of enviromental metabolites needed in new model
% blacklist_flag - zero for model reduction

% Usage 2: Blacklisting select proteins
% [model] = blacklist_tfs(fullmodel,blacklist)

% Input
% fullmodel - Model from which proteins are to be blacklisted
% blacklist - List of proteins to be blacklisted

% Usage 3: Arrange such that Protein = [protein;receptors]'
% [model] = blacklist_tfs(model,list,0);
 
% Input
% model - model in which arrangement is to be made
% list - empty cell array for rearrangement
% blacklist_flag - zero for rearrangement
%**************************************************************************
function [newmodel,nregp,nrecp] = blacklist_tfs(model,list,blacklist,recurmodel)
if nargin < 4
    newmodel.Gene = {};
    newmodel.Regulators = {};
    newmodel.Proteins = {};
    newmodel.Metabolites = {};
    newmodel.RS = sparse(0,0);
    newmodel.Coefficient = sparse(0,0);
    newmodel.trate = sparse(0,0);
    newmodel.brate = sparse(0,0);
    %newmodel.metabRS = sparse(0,0);
    newmodel.GeneRules = {};
    newmodel.mRNA = zeros(0,0);
    newmodel.Prot = zeros(0,0);
    newmodel.maxProt = zeros(0,0);
    newmodel.Metab = zeros(0,0);
    newmodel.maxMetab = zeros(0,0);
    newmodel.Reg = zeros(0,0);
    newmodel.maxReg = zeros(0,0);
else
    newmodel.Gene = recurmodel.Gene;
    newmodel.GeneRules = recurmodel.GeneRules;
    newmodel.Regulators = recurmodel.Regulators;
    newmodel.Proteins = recurmodel.Proteins;
    newmodel.Metabolites = recurmodel.Metabolites;
    newmodel.RS = recurmodel.RS;
    newmodel.Coefficient = recurmodel.Coefficient;
    newmodel.trate = recurmodel.trate;
    newmodel.brate = recurmodel.brate;
    %newmodel.metabRS = recurmodel.metabRS ;
    newmodel.mRNA = recurmodel.mRNA;
    newmodel.Prot = recurmodel.Prot;
    newmodel.maxProt = recurmodel.maxProt;
    newmodel.Metab = recurmodel.Metab;
    newmodel.maxMetab = recurmodel.maxMetab;
    newmodel.Reg = recurmodel.Reg;
    newmodel.maxReg = recurmodel.maxReg;
end
if nargin < 3
    blacklist = 1;
end
if isempty(list) && ~blacklist
    newmodel = model;
    %Arrange such that Protein = [protein;receptors]'
%     [recp,~] = find(model.metabRS);
%     [ngap,~] = find(model.trate);
%     newmodel.Gene = model.Gene;
%     newmodel.Protein = [model.Protein(ngap);model.Protein(recp)];
%     newmodel.Metabolite = model.Metabolite;
%     newmodel.GeneRules = model.GeneRules;
%     newmodel.brate = model.brate;
%     newmodel.mRNA = model.mRNA;
%     newmodel.RS = [model.RS(:,ngap),model.RS(:,recp)];
%     newmodel.Coefficient = [model.Coefficient(:,ngap),model.Coefficient(:,recp)];
%     newmodel.trate = [model.trate(ngap,:);model.trate(recp,:)];        
%     newmodel.metabRS = [model.metabRS(ngap,:);model.metabRS(recp,:)]; 
%     newmodel.Prot = [model.Prot(ngap);model.Prot(recp)];
%     newmodel.maxProt = [model.maxProt(ngap);model.maxProt(recp)];
%     nregp = length(ngap);
%     nrecp = length(recp);
    return
end

%Balcklist
if blacklist
    newmodel.GeneRules = model.GeneRules;
    newmodel.GeneRules(cellfun('isempty',newmodel.GeneRules)) = {'reg'};
    newmodel.Proteins = model.Proteins;
    newmodel.Gene = model.Gene;
    newmodel.Metabolites = model.Metabolites;
    newmodel.Regulators = model.Regulators;
    newmodel.RS = model.RS;
    newmodel.Coefficient = model.Coefficient;
    newmodel.trate = model.trate;
    newmodel.metabRS = model.metabRS;
    newmodel.brate = model.brate;   
end

listsize = length(list);
status = zeros(listsize,1);
 
remove_pind = zeros(length(model.Protein),1);
jlist = 1;
diag_prot = cell(listsize,1);
while jlist <= listsize  
    fprintf('\nAdditions Made to Model for Protein %s',list{jlist});
    mxind = strcmpi(model.Metabolite,list{jlist});
    if any(mxind)
        regpind = full(logical(model.metabRS(:,mxind)));
    else
        regpind = strcmpi(model.Protein,list{jlist});
    end
    if any(regpind)        
        gind = find(model.RS(:,regpind));
        if blacklist
            fprintf('\nBlacklisted Proteins: %s',model.Protein{regpind});            
            remove_pind(regpind) = 1;
        end
        kgene = 1;
        diag_gene = cell(length(gind),1);
        while kgene <= length(gind)
            run_flag = 0;
            if blacklist
                %Do blacklist duties          
                %Remove Regulator from all rules of affected Genes
                for jrule = 1:length(newmodel.GeneRules{gind(kgene)})
                    kind = strfind(newmodel.GeneRules{gind(kgene)}{jrule},list{jlist});
                    if ~isempty(kind)
                        GeneRules1 = newmodel.GeneRules{gind(kgene)}{jrule}(1:kind-1);
                        GeneRules2 = newmodel.GeneRules{gind(kgene)}{jrule}(kind+length(list{jlist})+1:end);
                        if ~isempty(GeneRules2)
                            newmodel.GeneRules{gind(kgene)}{jrule} = [GeneRules1 GeneRules2];                    
                        else
                            if jrule < length(newmodel.GeneRules{gind(kgene)})
                                newmodel.GeneRules{gind(kgene)}{jrule} = GeneRules1(1:end);
                            elseif jrule == length(newmodel.GeneRules{gind(kgene)})
                                newmodel.GeneRules{gind(kgene)}{jrule} = GeneRules1(1:end-1);
                            end
                        end
                    end
                end 
            else
                %generate_subsystem duties                   
                diag_gene{kgene} = model.Gene{gind(kgene)};
                diag_prot{jlist} = model.Protein{regpind};
                pind = logical(model.trate(:,gind(kgene)));
                new_pind = strcmpi(newmodel.Protein,model.Protein{pind});
                new_gind = strcmp(newmodel.Gene,model.Gene{gind(kgene)});%new Gene
                new_regpind = strcmpi(newmodel.Protein,model.Protein{regpind});%new Protein
                if any(mxind)
                    new_mxind = strcmp(newmodel.Metabolite,model.Metabolite{mxind});
                end
                if any(new_gind) && any(new_regpind)%=>Protein is there
                    %change values for
                    %a) regprotein
                    %b) gene protein
                    change_val(gind(kgene),pind,regpind,new_gind,new_pind,new_regpind);
                    newmodel.trate(new_regpind,new_gind) = model.trate(regpind,gind(kgene));
                    if any(mxind)
                        change_mval(regpind,mxind,new_regpind,new_mxind);                  
                    end 
                elseif any(new_gind)
                    %add new reg_protein
                    change_val(gind(kgene),pind,regpind,new_gind,new_pind,length(new_regpind)+1);
                    newmodel.Protein{length(new_regpind)+1} = model.Protein{regpind};  
                    newmodel.Prot(length(new_regpind)+1) = model.Prot(regpind);
                    newmodel.maxProt(length(new_regpind)+1) = model.maxProt(regpind);
                    newmodel.trate(length(new_regpind)+1,new_gind) = model.trate(regpind,gind(kgene));
                    if any(mxind)
                        change_mval(regpind,mxind,length(new_regpind)+1,new_mxind);
                    else
                        newmodel.metabRS(length(new_regpind)+1,:) =...
                        sparse(1,length(newmodel.Metabolite));
                    end                     
                elseif any(new_regpind)
                    change_val(gind(kgene),pind,regpind,length(new_gind)+1,length(new_pind)+1,new_regpind);
                    %add new gene
                    newmodel.Gene{length(new_gind)+1} = model.Gene{gind(kgene)};
                    newmodel.mRNA(length(new_gind)+1) = model.mRNA(gind(kgene));
                    %add its corresponding protein
                    newmodel.Protein{length(new_pind)+1} = model.Protein{pind};
                    newmodel.Prot(length(new_pind)+1) = model.Prot(pind);
                    newmodel.maxProt(length(new_pind)+1) = model.maxProt(pind);
                    newmodel.trate(length(new_pind)+1,length(new_gind)+1) = model.trate(pind,gind(kgene));
                    if any(mxind)
                        change_mval(regpind,mxind,new_regpind,new_mxind);              
                    end
                    newmodel.metabRS(length(new_pind)+1,:) = sparse(1,length(newmodel.Metabolite));
                else
                    change_val(gind(kgene),pind,regpind,length(new_gind)+1,length(new_pind)+1,length(new_regpind)+2);
                    %add new gene
                    newmodel.Gene{length(new_gind)+1} = model.Gene{gind(kgene)};
                    newmodel.mRNA(length(new_gind)+1) = model.mRNA(gind(kgene));%Concentration
                    %add it corresponding protein
                    newmodel.Protein{length(new_pind)+1} = model.Protein{pind};                    
                    newmodel.trate(length(new_pind)+1,length(new_gind)+1) = model.trate(pind,gind(kgene));
                    newmodel.metabRS(length(new_pind)+1,:) = sparse(1,length(newmodel.Metabolite));
                    newmodel.Prot(length(new_pind)+1) = model.Prot(pind);%Concentration
                    newmodel.maxProt(length(new_pind)+1) = model.maxProt(pind);
                    %add new reg_protein
                    newmodel.Protein{length(newmodel.Protein)+1} = model.Protein{regpind};
                    newmodel.Prot(length(newmodel.Protein)+1) = model.Prot(regpind);
                    newmodel.maxProt(length(newmodel.Protein)+1) = model.maxProt(regpind);
                    newmodel.trate(length(newmodel.Protein),length(new_gind)+1) = model.trate(regpind,gind(kgene));
                    if any(mxind)
                        change_mval(regpind,mxind,length(newmodel.Protein),new_mxind);
                    else
                        newmodel.metabRS(length(newmodel.Protein),:) = sparse(1,length(newmodel.Metabolite));
                    end 
                end                   
                
                if any(mxind)
                    more_regpind = logical(model.trate(:,gind(kgene)));
                    if any(regpind)                        
                        more_prot = model.Protein(more_regpind);                
                        if ~run_flag
                            recur.Gene = newmodel.Gene;
                            recur.GeneRules = newmodel.GeneRules;
                            recur.Protein = newmodel.Protein;
                            recur.Metabolite = newmodel.Metabolite;
                            recur.RS = newmodel.RS;
                            recur.Coefficient = newmodel.Coefficient;
                            recur.trate = newmodel.trate;
                            recur.brate = newmodel.brate;
                            recur.metabRS = newmodel.metabRS;
                            recur.mRNA = newmodel.mRNA;
                            recur.Prot = newmodel.Prot;
                            recur.maxProt = newmodel.maxProt;
                            [recur] = blacklist_tfs(model,more_prot,0,recur);
                            newmodel.Gene = recur.Gene;
                            newmodel.GeneRules = recur.GeneRules;
                            newmodel.Protein = recur.Protein ;
                            newmodel.Metabolite = recur.Metabolite;
                            newmodel.RS = recur.RS;
                            newmodel.Coefficient = recur.Coefficient;
                            newmodel.trate = recur.trate;
                            newmodel.brate = recur.brate;
                            newmodel.metabRS = recur.metabRS;
                            newmodel.mRNA = recur.mRNA;
                            newmodel.Prot = recur.Prot;
                            newmodel.maxProt = recur.maxProt;
                        end
                    end
                end                 
            end
            kgene = kgene + 1;
        end  
    else
        fprintf('\nProtein %s Not Present in model',list{jlist});        
    end    
    jlist = jlist + 1;
end
%Print Diagnostics



%Blaclist
remove_pind = logical(remove_pind);
if blacklist
    if any(remove_pind)
        newmodel.RS(:,remove_pind) = [];
        newmodel.Coefficient(:,remove_pind) = [];
        newmodel.metabRS(remove_pind,:) = [];
        newmodel.trate(remove_pind,:) = [];
        newmodel.Prot(remove_pind) = [];
        newmodel.maxProt(remove_pind) = [];
        nonzero_pind = find(remove_pind);
        for iremove = 1:length(nonzero_pind)
            newmodel.Protein{nonzero_pind(iremove)} = {};
        end
        newmodel.Protein = newmodel.Protein(~cellfun('isempty',newmodel.Protein));
    end 
else
    %Generate Subsystem
    fprintf('\n-------Additions Made to Model for Protein %s-------',list{1:listsize});
    fprintf('\n\t\t %s \t\t\t %s \t\t\t %s',diag_gene{1:length(diag_gene)});
    newmodel.Gene = newmodel.Gene';
    newmodel.GeneRules = newmodel.GeneRules';
    newmodel.Protein = newmodel.Protein';
    newmodel.Metabolite = newmodel.Metabolite';
    newmodel.RS = newmodel.RS;
    newmodel.Coefficient = newmodel.Coefficient;
    newmodel.trate = newmodel.trate;
    newmodel.brate = newmodel.brate';
    newmodel.metabRS = newmodel.metabRS;
    newmodel.Prot = newmodel.Prot';
    newmodel.maxProt = newmodel.maxProt';
    newmodel.mRNA = newmodel.mRNA';
    nregp = length(find(newmodel.metabRS));
    nrecp = length(find(newmodel.trate));
end

function change_mval(pind,mxind,new_pind,new_mxind)
    if any(new_mxind)
        newmodel.metabRS(new_pind,new_mxind) = model.metabRS(pind,mxind);
    else
        newmodel.metabRS(new_pind,length(new_mxind)+1) = model.metabRS(pind,mxind);
        newmodel.Metabolite{length(new_mxind)+1} = model.Metabolite{mxind};        
    end 
end
function change_val(old_gindx,old_pindx,old_regpindx,new_gindx,new_pindx,new_regpindx)
    %RS    
    newmodel.RS(new_gindx,new_pindx) = model.RS(old_gindx,old_pindx);
    newmodel.RS(new_gindx,new_regpindx) = model.RS(old_gindx,old_regpindx);
    %Coefficient
    newmodel.Coefficient(new_gindx,new_pindx) = model.Coefficient(old_gindx,old_pindx);
    newmodel.Coefficient(new_gindx,new_regpindx) = model.Coefficient(old_gindx,old_regpindx);      
    %brate
    newmodel.brate(new_gindx) = model.brate(old_gindx); 
    %Concentrations
    newmodel.mRNA(new_gindx) = model.mRNA(old_gindx);
    newmodel.Prot(new_pindx) = model.Prot(old_pindx);
    newmodel.Prot(new_regpindx) = model.Prot(old_regpindx);
    newmodel.maxProt(new_pindx) = model.maxProt(old_pindx);
    newmodel.maxProt(new_regpindx) = model.maxProt(old_regpindx);
    %Regulatory Rules
    if ~islogical(new_gindx)
        if length(newmodel.GeneRules) > new_gindx
            if ~isempty(newmodel.GeneRules{new_gindx})   
                if ~iscell(newmodel.GeneRules{new_gindx})
                    newmodel.GeneRules{new_gindx} = changeGeneRules(newmodel.GeneRules(new_gindx),model.Protein{old_regpindx});
                else
                    newmodel.GeneRules{new_gindx} = changeGeneRules(newmodel.GeneRules{new_gindx},model.Protein{old_regpindx});
                end
            else
                newmodel.GeneRules{new_gindx} = model.Protein{old_regpindx};
            end
        else
            newmodel.GeneRules{new_gindx} = model.Protein{old_regpindx};
        end
    else
        if ~isempty(newmodel.GeneRules{new_gindx})
            if ~iscell(newmodel.GeneRules{new_gindx})
                newmodel.GeneRules{new_gindx} = changeGeneRules(newmodel.GeneRules(new_gindx),model.Protein{old_regpindx});
            else
                newmodel.GeneRules{new_gindx} = changeGeneRules(newmodel.GeneRules{new_gindx},model.Protein{old_regpindx});
            end
        else
            newmodel.GeneRules{new_gindx} = model.Protein{old_regpindx};
        end
    end
end

end