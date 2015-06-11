%function newmodel = separatesubsystems(trnmodel,regprotein)
%creating a subsytem of genes and/or TFs for evaluation disregarding the
%entire system/network
function [newmodel,newregprot] = separatesubsystems(gene,model,regprotein)
% ngene = length(model.Gene);
% tnreg = length(model.Protein);
% [coefficient,srate] = parameter_return(model.allpar,model,ngene,tnreg);
coefficient = model.Coefficient;
srate = model.srate;

Protein = model.Protein;
Metabolite = model.Metabolite;

%gene = {'aceA'};
newregprot = {};
rgene = {};
rprotein = {};


for ig = 1:length(gene)    
    genetf = strcmp(gene{ig},model.Gene);
    recursive_func(genetf);     
end
gene = [gene;rgene];

if ~isempty(gene)
    for iprot = 1:length(rprotein)
        prot_tf = strcmp(rprotein{iprot},regprotein);
        if any(prot_tf)
            newregprot = [newregprot;rprotein{iprot}];
        end
    end

    newg_count = length(gene);
    newp_count = length(rprotein);

    newmodel = struct();
    newmodel.Gene = cell(newg_count,1);
    newmodel.GeneRules = cell(newg_count,1);
    newmodel.Protein = cell(newp_count,1);
    newmodel.Metabolite = model.Metabolite;
    newmodel.RS = sparse(newg_count,newp_count);
    newmodel.brate = zeros(newg_count,1);
    newmodel.Coefficient = sparse(newg_count,newp_count);
    newmodel.trate = sparse(length(newregprot),newg_count);
    newmodel.srate = sparse(newg_count,newp_count);
    newmodel.Kmax = sparse(length(rprotein),length(Metabolite));


    for ig = 1:length(gene)
        gtf = strcmp(gene{ig},model.Gene);
        if any(gtf)
            newmodel.Gene{ig} = gene{ig};
            newmodel.GeneRules{ig} = model.GeneRules{logical(gtf)}; 
            newmodel.brate(ig) = model.brate(gtf);
            for ip = 1:length(rprotein)
                ptf = strcmp(rprotein{ip},model.Protein);
                if any(ptf)
                    newmodel.Protein{ip} = rprotein{ip};
                    if logical(model.RS(gtf,ptf))
                        newmodel.RS(ig,ip) = model.RS(gtf,ptf);
                    end
                    if logical(coefficient(gtf,ptf))
                        newmodel.Coefficient(ig,ip) = coefficient(gtf,ptf);
                    end
                    if logical(srate(gtf,ptf))
                        newmodel.srate(ig,ip) = srate(gtf,ptf);
                    end
                    prtf = strcmp(rprotein{ip},regprotein);
                    if logical(model.trate(prtf,gtf))
                        new_prtf = strcmp(rprotein{ip},newregprot);
                        newmodel.trate(new_prtf,ig) = model.trate(prtf,gtf);
                    end
                    newmodel.Kmax(ip,:) = model.Kmax(ptf,:);                
                end
            end
        end
    end     
end

newmodel.allpar = parameter_vector(newmodel.Coefficient,newmodel.srate,newg_count);
newmodel = rmfield(newmodel,{'Coefficient','srate'});



function recursive_func(indx)                    
    regindx = model.RS(logical(indx),:);
    reg = model.Protein(logical(regindx));
    for ireg = 1:length(reg)
        regprot_tf = strcmp(reg{ireg},regprotein);
        if any(regprot_tf)
            prot_gindx = model.trate(regprot_tf,:);
            rgtf = strcmp(model.Gene{logical(prot_gindx)},rgene);
            if ~any(rgtf)
                rgene = [rgene;model.Gene{logical(prot_gindx)}];
                recursive_func(prot_gindx);
            end
        else            
            metabprot_tf = strcmp(reg{ireg},Protein);
            if any(metabprot_tf)
                if ~any(strcmp(Protein{logical(metabprot_tf)},rprotein))
                    rprotein = [rprotein;Protein{logical(metabprot_tf)}];
                end
            end
        end
        if ~any(strcmp(reg{ireg},rprotein));
            rprotein = [rprotein;reg{ireg}];
        end
    end       
end
end


       
        
