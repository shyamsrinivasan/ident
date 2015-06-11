%Function to obtain smaller self contained sub-systems from the complete
%trnmodel
%Specify only the regulating metabolites that are needed.
function [newmodel,newregprot,ngene,tnreg,ngap] = generate_subsystem(model,metab,regprotein)
% metab = {'Oxygen'};
newGene = {};
protein = {};
newRS = sparse(0,0);
newCoeff = [];
%newsrate = [];
newtrate = [];
newmetabRS = [];


for im = 1:length(metab)
    %search for corresponding receptor 
    mind = strcmp(metab{im},model.Metabolite);
    [recp,~] = find(model.metabRS(:,mind));        
    kprot = length(protein);
    
    irecp = 1;
    while irecp <= length(recp)
        %if kprot == 0             
            kprot = kprot + 1;                  
        %end
        protein{kprot} = model.Protein{recp(irecp)};        
        newmetabRS(kprot,im) = 1;             
        [reg_gene,~] = find(model.RS(:,recp(irecp)));  
        
        ireg_gene = 1;
        while ireg_gene <= length(reg_gene)            
            %identify protein of regulatory gene
            reg_prot = find(model.trate(:,reg_gene(ireg_gene)));            
            for ireg = 1:length(reg_prot)
                newprot_indx = strcmp(model.Protein{reg_prot(ireg)},protein);
                if ~any(newprot_indx)%protein not in new system
                    kprot = kprot + 1;
                    protein{kprot} = model.Protein{reg_prot(ireg)};
                    newmetabRS(kprot,im) = 0;    
                    
                    %determine genes regulated by new protein
                    [gind,~] = find(model.RS(:,reg_prot(ireg)));
                    for kgind = 1:length(gind)
                        new_gind = strcmp(model.Gene{gind(kgind)},newGene);
                        kg = length(newGene);
                        if ~any(new_gind)%gene not in newGene => add it
                            newGene{kg+1} = model.Gene{gind(kgind)};
                            change_val(gind(kgind),reg_prot(ireg),...
                                       kg+1,kprot)
                            newtrate(kprot,kg+1) =...
                            model.trate(reg_prot(ireg),gind(kgind));
                        else%gene in newGene
                            change_val(gind(kgind),reg_prot(ireg),...
                                       new_gind,kprot)%                             
                            newtrate(kprot,new_gind) =...
                            model.trate(reg_prot(ireg),gind(kgind));%                            
                        end                            
                    end                 
                    
                    gene_indx = strcmp(model.Gene{reg_gene(ireg_gene)},newGene);
                    prot_indx = strcmp(model.Protein{recp(irecp)},protein);
                    kg = length(newGene);
                    if any(gene_indx)
                        change_val(reg_gene(ireg_gene),recp(irecp),...
                                   gene_indx,prot_indx)%                         
                        newtrate(kprot,gene_indx) =...
                        model.trate(reg_prot(ireg),reg_gene(ireg_gene));                       
                        
                    else
                        newGene{kg+1} = model.Gene{reg_gene(ireg_gene)};
                        change_val(reg_gene(ireg_gene),recp(irecp),...
                                   kg+1,prot_indx)%                         
                        newtrate(kprot,kg+1) = model.trate(reg_prot(ireg),reg_gene(ireg_gene));
                        
                    end
                else%protein already in new system
                    fprintf('place holder');
                    %determine genes regulated by new protein
                    [gind,~] = find(model.RS(:,reg_prot(ireg)));
                    for kgind = 1:length(gind)
                        new_gind = strcmp(model.Gene{gind(kgind)},newGene);
                        kg = length(newGene);
                        if ~any(new_gind)%gene not in newGene => add it
                            newGene{kg+1} = model.Gene{gind(kgind)};
                            change_val(gind(kgind),reg_prot(ireg),...
                                       kg+1,newprot_indx)
                            newtrate(newprot_indx,kg+1) =...
                            model.trate(reg_prot(ireg),gind(kgind));
                        else%gene in newGene
                            change_val(gind(kgind),reg_prot(ireg),...
                                       new_gind,newprot_indx)%                             
                            newtrate(newprot_indx,new_gind) =...
                            model.trate(reg_prot(ireg),gind(kgind));%                            
                        end
                    end
                    
                    gene_indx = strcmp(model.Gene{reg_gene(ireg_gene)},newGene);
                    prot_indx = strcmp(model.Protein{recp(irecp)},protein);
                    kg = length(newGene);
                    if any(gene_indx)
                        change_val(reg_gene(ireg_gene),recp(irecp),...
                                   gene_indx,prot_indx)%                         
                        newtrate(newprot_indx,gene_indx) =...
                        model.trate(reg_prot(ireg),reg_gene(ireg_gene));                       
                        
                    else
                        newGene{kg+1} = model.Gene{reg_gene(ireg_gene)};
                        change_val(reg_gene(ireg_gene),recp(irecp),...
                                   kg+1,prot_indx)%                         
                        newtrate(newprot_indx,kg+1) = model.trate(reg_prot(ireg),reg_gene(ireg_gene));
                        
                    end
                    
                end
                
                %adding of additional regulators follows here
                kprot = length(protein);
                for ig = 1:length(newGene)
                    gind_1 = strcmp(newGene{ig},model.Gene);
                    [prot,~] = find(model.trate(:,gind_1));
                    iprot = 1;
                    while iprot <= length(prot)  
                        [generow,~] = find(model.RS(:,prot(iprot)));
                        if ~any(strcmp(model.Protein{prot(iprot)},protein))
                            %Add protein to the list and the model                            
                            protein{kprot+1} = model.Protein{prot(iprot)};                         
                            %identify genes in newGene regulated by new protein                     
                            for ign = 1:length(generow)
                                gcomp = strcmp(model.Gene{generow(ign)},newGene);                
                                if any(gcomp)%  
                                    change_val(generow(ign),prot(iprot),...
                                               gcomp,kprot+1);
                                    newtrate(kprot+1,gcomp) =...
                                    model.trate(prot(iprot),generow(ign));%                                                  
                                    newmetabRS(kprot+1,im) = 0;
                                    
                                end
                            end 
                            kprot = kprot + 1;
                        else%else case is not required since proteins that re there are
                            %already accounted for 
                            
                            %change existing matrices            
                %             pind = strcmp(trnmodel.Protein{prot(iprot)},protein);
                % %             comp_prot(pind);
                %             %[generow,~] = find(trnmodel.RS(:,prot(iprot)));
                %             for ign = 1:length(generow)
                %                 gcomp = strcmp(trnmodel.Gene{generow(ign)},newGene);
                %                 if any(gcomp)
                %                     newRS(gcomp,pind) = trnmodel.RS(generow(ign),prot(iprot));
                %                     newCoeff(gcomp,pind) = trnmodel.Coefficient(generow(ign),prot(iprot));
                %                     newsrate(gcomp,pind) = trnmodel.srate(generow(ign),prot(iprot));
                %                     newtrate(pind,gcomp) = trnmodel.trate(prot(iprot),generow(ign));
                %                 end
                %             end            

                        end                
                        iprot = iprot + 1;
                    end            
                end                 
            end            
            ireg_gene = ireg_gene + 1;
        end
        irecp = irecp + 1;       
    end
end

newRS = sparse(newRS);
newCoeff = sparse(newCoeff);
%newsrate = sparse(newsrate);
newtrate = sparse(newtrate);
protein = protein';
newmetabRS = sparse(newmetabRS);

ngene = length(newGene);
tnreg = length(protein);
newregprot = protein(cellfun('isempty',regexp(protein,'\w+.?Recp')));
ngap = length(newregprot);

newmodel.Gene = newGene';
newmodel.Protein = cell(ngap,1);
newmodel.Metabolite = metab;
newmodel.GeneRules = cell(ngene,1);


newmodel.RS = sparse(ngene,ngap);
newmodel.Coefficient = sparse(ngene,ngap);
%newmodel.srate = sparse(ngene,ngap);
newmodel.trate = sparse(ngap,ngene);
newmodel.brate = newbrate;
newmodel.metabRS = sparse(ngap,length(metab));


reg_diff = setdiff(model.Protein,protein);
model1 = blacklist_tfs(model,reg_diff,regprotein);

for ig = 1:length(newmodel.Gene)
    gene_ind = strcmp(newmodel.Gene{ig},model.Gene);
    if any(gene_ind)
        newmodel.GeneRules{ig,1} = model1.GeneRules{gene_ind};
    else
        newmodel.GeneRules{ig,1} = {};
    end
end


jprot = 1;
indx = zeros(1,tnreg);
for igap = 1:ngap
    gapindx = strcmp(newregprot{igap},protein);
    newmodel.Protein{jprot,1} = protein{gapindx};
    newmodel.RS(:,jprot) = newRS(:,gapindx);
    newmodel.Coefficient(:,jprot) = newCoeff(:,gapindx);
    %newmodel.srate(:,jprot) = newsrate(:,gapindx);
    newmodel.trate(jprot,:) = newtrate(gapindx,:);    
    newmodel.metabRS(jprot,:) = newmetabRS(gapindx,:);
    
    indx(gapindx) = 1;  
    jprot = jprot + 1;
end

newmodel.Protein = [newmodel.Protein;protein(~indx)];
newmodel.RS = [newmodel.RS,newRS(:,~indx)];
newmodel.Coefficient = [newmodel.Coefficient,newCoeff(:,~indx)];
%newmodel.srate = [newmodel.srate,newsrate(:,~indx)];
newmodel.metabRS = [newmodel.metabRS; newmetabRS(~indx,:)];

newmodel.allpar = parameter_vector(newmodel.Coefficient,newmodel.srate,ngene);
fields = {'Coefficient','srate'};
newmodel = rmfield(newmodel,fields); 


function change_val(old_geneindx,old_protindx,new_geneindx,new_protindx)
    newRS(new_geneindx,new_protindx) = model.RS(old_geneindx,old_protindx);
    newCoeff(new_geneindx,new_protindx) = model.Coefficient(old_geneindx,old_protindx);
    %newsrate(new_geneindx,new_protindx) = model.srate(old_geneindx,old_protindx);    
    newbrate(new_geneindx) = model.brate(old_geneindx);    
end



end