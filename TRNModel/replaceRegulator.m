% function [model] = replaceRegulator(model,oldreg,newreg)
% Substitute transcriptional regulator oldreg with newreg in all instances
function [model] = replaceRegulator(model,oldreg,newreg)
tfr_old = strcmpi(oldreg,model.Regulators);
%Find where oldreg is being used
gind = find(model.RS(:,tfr_old)~=0);
for ig = 1:length(gind)
    %Split Rule Terms
    nrules = length(model.GeneRules{gind(ig)});
    if nrules == 1
        try
            [terms,~,andpos,orterms,andterms] =...
            return_terms(model.GeneRules{gind(ig)});
        catch
            [terms,~,andpos,orterms,andterms] =...
            return_terms(model.GeneRules{gind(ig)}{1});
        end
        %Changes to Gene Rules
        model.GeneRules{gind(ig)} =...
        adjustRules(oldreg,newreg,terms,orterms,andterms,andpos);       
        %Changes to Other Elements of model
        model = adjustmodel(model,oldreg,newreg);        
    elseif nrules > 1
        irule = 1;
        while irule <= nrules 
            [terms,~,andpos,orterms,andterms] =...
            return_terms(model.GeneRules{gind(ig)}{irule}); 
            %Changes to Gene Rules
            model.GeneRules{gind(ig)}{irule} =...
            adjustRules(oldreg,newreg,terms,orterms,andterms,andpos);       
            %Changes to Other Elements of model
            model = adjustmodel(model,oldreg,newreg);   
            irule = irule + 1;
        end
    end   
end
return 

function newRule = adjustRules(oldreg,newreg,terms,orterms,andterms,andpos)
    old_pr = strcmpi(oldreg,terms);
    new_pr = strcmpi(newreg,terms);
    %If oldreg is present as a regulating element
    if any(old_pr) && ~any(new_pr)
        %Change oldreg to newreg
        terms{old_pr} = newreg;
        %change Gene Rules
        %Determine connection logic
        if ~isemptyr(orterms) && isempty(andpos)
            %OR '|'
            newRule = strjoin(terms,'|');
        elseif ~isemptyr(andterms) 
            %AND '&'
            newRule = strjoin(terms,'&');
        end            
    end
return

function [model] = adjustmodel(model,oldreg,newreg)
    kprot = length(model.Enzyme)+1;
    jmetab = length(model.Metabolites)+1;
    jreg = length(model.Regulators)+1;

    old_tfp = strcmpi(oldreg,model.Enzyme);
    old_tfm = strcmpi(oldreg,model.Metabolites);
    old_tfr = strcmpi(oldreg,model.Regulators);
    new_tfp = strcmpi(newreg,model.Enzyme);
    new_tfm = strcmpi(newreg,model.Metabolites);
    new_tfr = strcmpi(newreg,model.Regulators);
    if any(old_tfp) && ~any(new_tfp)
        model.Enzyme{kprot} = model.Enzyme{old_tfp};
        model.Enzyme{old_tfp} = newreg;  
        model.SSprot(kprot) = model.SSprot(old_tfp);
        model.pmRatio(kprot) = model.pmRatio(old_tfp);
        model.Vss(kprot) = model.Vss(old_tfp);
        model.Kcat(kprot) = model.Kcat(old_tfp);
        model.delSGr(kprot) = model.delSGr(old_tfp);
        model.Keq(kprot) = model.Keq(old_tfp);
        model.S(:,kprot) = 0;
        model.SI(:,kprot) = 0;
        model.pmeter.K(:,kprot) = 0;
        model.pmeter.KI(:,kprot) = 0;
        model.pmeter.KIact(:,kprot) = 0;
        model.pmeter.KIihb(:,kprot) = 0;
%         kprot = kprot + 1;
    end
    if any(old_tfm) && ~any(new_tfm)
        model.Metabolites{jmetab} = model.Metabolites{old_tfm};
        model.Metabolites{old_tfm} = newreg;
        model.Metab(jmetab) = 0;
        model.maxMetab(jmetab) = 0;
        model.S(jmetab,:) = 0;
        model.pmeter.K(jmetab,:) = 0;
        model.SI(jmetab,:) = 0;
        model.pmeter.KIact(jmetab,:) = 0;
        model.pmeter.KIihb(jmetab,:) = 0;
        model.MC(jmetab) = model.MC(old_tfm);
        model.MClow(jmetab) = model.MClow(old_tfm);
        model.MChigh(jmetab) = model.MChigh(old_tfm);
%         jmetab = jmetab + 1;
    end
    if any(old_tfr) && ~any(new_tfr)
        model.Regulators{jreg} = model.Regulators{old_tfr};
        model.Regulators{old_tfr} = newreg;
        model.RS(:,jreg) = 0;
        model.Coefficient(:,jreg) = 0; 
        model.trate(jreg,:) = 0;
        model.Kb(jreg) = 0;
        model.Kub(jreg) = 0;
%         model.SSreg(jreg) = 0;
        model.relReg(jreg) = 0;
%         jreg = jreg + 1;
    end
return

function [terms,orpos,andpos,orterms,andterms] = return_terms(GeneRules)
%Determine all the proteins in the rule
    orpos = strfind(GeneRules,'|');
    andpos = strfind(GeneRules,'&');  
    orterms = strsplit(GeneRules,'|');
    andterms = strsplit(GeneRules,'&');

    if ~isempty(orpos) && isempty(andpos)
        terms = strsplit(GeneRules,'|');
    elseif ~isempty(andpos) && isempty(orpos)
        terms = strsplit(GeneRules,'&');
    elseif isempty(orpos) && isempty(andpos)
        terms = {GeneRules};
    end    

% 
% nprots = length(terms);  
% new_name = cell(nprots,1);
% new_baff = zeros(nprots,1);
% new_sign = zeros(nprots,1);
% for iprot = 1:nprots
%     tfp = strcmp(terms{iprot},prot_name);
%     if any(tfp)
%         new_name{iprot} = prot_name{tfp};
%         new_baff(iprot) = bindaff(tfp);
%         new_sign(iprot) = RS(1,tfp);                    
%     end
% end
% inh_ind = logical(new_sign < 0);%inhibiting indices
% act_ind = logical(new_sign > 0);%activating indices

return