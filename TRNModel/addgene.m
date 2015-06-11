%function [trnmodel,gindx,T] = addgene(trnmodel,regprotein,rgene,rmetab,...
%                                      rprot,defparval,basis_flag)
%July 2014 version 2.0
%**************************************************************************
% ************Input
% trnmodel -    TRN model to which additional genes and proteins are added
% regprotein -  List of proteins that behave as TFs
% rgene -       Genes corresponding to proteins in regprotein
% rmetab -      Regulatory interactions that are observed as a result of
%               metabolite (not a TF)
% rprot -       Regulatory interactions that are observed as a result of a
%               protein (a TF)
% defparval -   A structure of default parameter values
%***********Optional Input
% basis_flag - 
% **********Output
% trnmodel - TRN model with added genes, metabolites and proteins and their
%           regulatory interactions
%**************************************************************************
function [model] = addgene(model,D,defparval,basis_flag)
                                  
if nargin < 7
    basis_flag = true;
end
% 
rprot = D{4};
regprotein = D{1};
rgene = D{2};%genes corresponding to regulatory proteins  
rmetab = D{3};
nregprotein = length(D{1});
nrgene = length(rgene);
model.Metabolite = {};
%model.trate = sparse(nregprotein,length(model.Gene));

if ~isemptyr(rmetab)
    model.metabRS = sparse(0,0);              
end

deftrate = defparval.trate; %4.94e-6;
defbrate = defparval.brate;%1.66e-7;
defkmax = defparval.kmax;%0.5;
defks = defparval.ks;%0.1;

if nrgene == nregprotein    
    irm = length(model.Metabolite) + 1; 
    jrp = length(model.Protein) + 1;
    jrg = length(model.Gene) + 1;
    
    %Rates & Coefficients for added genes
    for irgene = 1:nrgene               
        rpindx = strcmpi(regprotein{irgene},model.Protein);
        rgene_mul = strrep(strsplit(rgene{irgene},','),'"','');
        for imul = 1:length(rgene_mul)            
            rgindx = strcmp(rgene_mul{imul},model.Gene);
            if any(rpindx) && any(rgindx)
                model.trate(rpindx,rgindx) = deftrate;
                model.Prot(rpindx) = D{6}(irgene);
                model.maxProt(rpindx) = D{7}(irgene);
            elseif any(rpindx)
                model.trate(rpindx,jrg) = deftrate;
                model.Gene{jrg} = rgene_mul{imul};
                model.GeneRules{jrg} = {};
                model.Regulator{jrg} = {};
                model.RS = [model.RS;sparse(1,length(rpindx))];
                model.Coefficient = [model.Coefficient;sparse(1,length(rpindx))];
                model.brate(jrg) = defbrate;
                model.mRNA(jrg) = D{5}(irgene);%mRNA concentration
                jrg = jrg + 1;
            elseif any(rgindx)
                model.trate(jrp,rgindx) = deftrate;
                model.Protein{jrp} = regprotein{irgene};
                model.RS = [model.RS,sparse(length(rgindx),1)];
                model.Coefficient = [model.Coefficient,sparse(length(rgindx),1)];
                model.Prot(jrp) = D{6}(irgene);%Relative Protein Concentration
                model.maxProt(jrp) = D{7}(irgene);%Maximum Protein Concentration
                jrp = jrp + 1;
            else
                model.trate(jrp,jrg) = deftrate;
                model.Protein{jrp} = regprotein{irgene};
                model.Gene{jrg} = rgene_mul{imul};
                model.GeneRules{jrg} = {};
                model.Regulator{jrg} = {};
                model.RS = [model.RS,sparse(length(rgindx),1)];
                model.RS = [model.RS;sparse(1,length(rpindx)+1)];
                model.Coefficient = [model.Coefficient,sparse(length(rgindx),1)];
                model.Coefficient = [model.Coefficient;sparse(1,length(rpindx)+1)];
                model.brate(jrg) = defbrate; 
                model.mRNA(jrg) = D{5}(irgene);
                model.Prot(jrp) = D{6}(irgene);
                model.maxProt(jrp) = D{7}(irgene);
                jrp = jrp + 1;
                jrg = jrg + 1;
            end
    
            %Gene-Metabolite regulatory Info
            if ~isempty(rmetab{irgene})
                ormetab = strsplit(rmetab{irgene},'|');
                andmetab = strsplit(rmetab{irgene},'&');
                if ~isempty(ormetab)
                    [model,jrg,jrp,irm] = addterms(model,ormetab,jrg,jrp,irm,'rmetab');
                elseif ~isempty(andmetab)
                    [model,jrg,jrp,irm] = addterms(model,andmetab,jrg,jrp,irm,'rmetab');
                end
            end
            
            %Gene-Protein regulatory Info
            if ~isempty(rprot{irgene})
                orprot = strsplit(rprot{irgene},'|');
                andprot = strsplit(rprot{irgene},'&');
                if ~isemptyr(orprot)
                    [model,jrg,jrp,irm] = addterms(model,orprot,jrg,jrp,irm,'rprot');
                elseif ~isemptyr(andprot)
                    [model,jrg,jrp,irm] = addterms(model,andprot,jrg,jrp,irm,'rprot');
                end 
            end
        end
    end  
end
model.Metabolite = (model.Metabolite)';

function [model,krg,krp,irm] = addterms(model,terms,krg,krp,irm,call_sign)
    jterm = 1;
    while jterm <= length(terms)
        [match] = regexp(terms{jterm},'(\w+.?)\[(\W+.?)\]','tokens');
        if isempty(match)
            return
        end
        if ~isempty(rgene{irgene})
            gindx = strcmp(rgene{irgene},model.Gene);
        end
        rmindx = strcmp(match{1}{1},model.Metabolite);
        rpindx = strcmpi(match{1}{1},model.Protein);
        switch match{1}{2}
            case '+'
                stoich = 1;
                coeff = defparval.accoeff;
            case '-'
                stoich = -1;
                coeff = defparval.repcoeff;
        end
        if (any(gindx) && any(rmindx)) || (any(gindx) && any(rpindx))            
            if strcmp(call_sign,'rmetab')
                pmindx = strcmp(sprintf('%sRecp',model.Metabolite{rmindx}),...
                            model.Protein); 
                newindx = pmindx;                
                newreg = sprintf('%sRecp',match{1}{1});
                model.metabRS(pmindx,rmindx) = 1;
            elseif strcmp(call_sign,'rprot')
                newindx = rpindx;
                newreg = match{1}{1};                
            end                                   

            if model.RS(gindx,newindx) == 0
                if ~isempty(model.GeneRules{gindx})
                    model.GeneRules{gindx} =...
                    changeGeneRules(model.GeneRules{gindx},...
                    sprintf('%s',model.Protein{newindx}));                                 
                else                                
                    model.GeneRules{gindx} = {newreg};
                end
                model.RS(gindx,newindx) = stoich; 
                model.Regulator{gindx}{end+1} = model.Protein{newindx}; 
                model.Coefficient(gindx,newindx) = coeff;
            end         
        elseif any(gindx)
            if strcmp(call_sign,'rmetab')
                model.Protein{krp} = sprintf('%sRecp',match{1}{1});
                model.Metabolite{irm} = match{1}{1};          
                model.metabRS(krp,irm) = 1;  
                irm = irm + 1;
            elseif strcmp(call_sign,'rprot')              
                model.Protein{krp} = match{1}{1};
            end                                                                      
            if ~isempty(model.GeneRules{gindx})
                model.GeneRules{gindx} =...
                changeGeneRules(model.GeneRules{gindx},model.Protein{krp});                             
            else
                model.GeneRules{gindx} = model.Protein{krp};
            end
            model.trate = [model.trate;sparse(1,length(gindx))];
            model.RS(gindx,krp) = stoich;
            model.Regulator{gindx}{end+1} = model.Protein{krp};
            model.Coefficient = [model.Coefficient,...
                                 sparse(find(gindx),1,coeff,length(gindx),1)]; 
            model.Prot(krp) = 0;%Protein concentration
            model.maxProt(krp) = 0;
            krp = krp + 1;                    
        elseif any(rmindx) || any(rpindx)
            if strcmp(call_sign,'rmetab')
                pmindx = strcmp(sprintf('%sRecp',model.Metabolite{rmindx}),...
                         model.Protein);
                newindx = pmindx;       
                model.metabRS(pmindx,rmindx) = 1; 
            elseif strcmp(call_sign,'rprot')
                newindx = rpindx;                
            end              
            model.Gene{krg} = rgene{irgene};
            model.GeneRules{krg} = model.Protein{newindx};  
            model.RS(krg,newindx) = stoich;  
            model.Regulator{krg} = model.Protein{newindx};
            model.Coefficient = [model.Coefficient;...
                                 sparse(1,find(newindx),coeff,1,length(newindx))];
            model.mRNA(krg) = 0;%mRNA concentration
            krg = krg + 1;       
        else
            if strcmp(call_sign,'rmetab')
                model.Protein{krp} = sprintf('%sRecp',match{1}{1});
                model.Metabolite{irm} = match{1}{1};
                model.metabRS(krp,irm) = 1;
                irm = irm + 1;
            elseif strcmp(call_sign,'rprot')
                model.Protein{krp} = match{1}{1};
            end         
            model.Gene{krg} = rgene{irgene};
            model.GeneRules{krg} =  model.Protein{krp}; 
            model.RS(krg,krp) = stoich;  
            model.Regulator{krg} =  model.Protein{krp};
            model.Coefficient = [model.Coefficient,...
                                 sparse(length(gindx),1)];
            model.Coefficient = [model.Coefficient;...
                                 sparse(1,krp,coeff,1,krg)];
            model.mRNA(krg) = 0;
            model.Prot(krp) = 0;
            model.maxProt(krp) = 0;
            krg = krg + 1;
            krp = krp + 1;            
        end
        jterm = jterm + 1;
    end 
end


end