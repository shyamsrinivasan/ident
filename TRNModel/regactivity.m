% function [pactnr,pactdr] = regactivity(prot_name,RS,Rules,bindaff,data)
% Function to calculate the regulatory activity of a promoter given binding
% affinities of trasncription factors
% Binding affinity = [TF]/K where K is the binding constant for activation
% or inhibition
function [pact_nr,pact_dr] = regactivity(prot_name,RS,Rules,bindaff,data)

pact = zeros(1,1);
if iscell(Rules)
    nrules = length(Rules);
else
    nrules = 1;
end
connectlogic = cell(nrules,1);
pactnr = zeros(nrules,1);
pactdr = zeros(nrules,1);

if nrules == 1
    if ~isempty(Rules)
        if iscell(Rules)            
            [terms,nprots,new_baff,inh_ind,act_ind] = returnterms_activity(Rules{1});
        else
            [terms,nprots,new_baff,inh_ind,act_ind] = returnterms_activity(Rules);
        end
        terms = terms(~cellfun('isempty',terms));
    %   Determine pactivitiy
        if iscell(Rules)
            [pactnr,pactdr] = pregactivity(new_baff,act_ind,inh_ind,Rules{1},data);
        else
            [pactnr,pactdr] = pregactivity(new_baff,act_ind,inh_ind,Rules,data);
        end             
    end
elseif nrules > 1
    irule = 1;
    while irule <= nrules  
        [terms,nprots,new_baff,inh_ind,act_ind] = returnterms_activity(Rules{irule});
        %Determine connectlogic between rules
        if nprots == length(orpos) + length(andpos);
            if ~isempty(orpos) && isemptyr(connectlogic{irule})
                connectlogic{irule} = Rules{irule}(orpos(end));
                Rules{irule} = Rules{irule}(1:end-1);
                terms{end}{2} = {};
            elseif ~isempty(andpos) && isemptyr(connectlogic{irule})
                connectlogic{irule} = Rules{irule}(andpos(end));
                Rules{irule} = Rules{irule}(1:end-1);
                terms{end}{2} = {};
            else
                connectlogic{irule} = {};
            end        
        end
        terms{end} = terms{end}(~cellfun('isempty',terms{end}));       
        %Determine pactivity
        [pactnr(irule),pactdr(irule)] = pregactivity(new_baff,act_ind,inh_ind,Rules{irule},data);
        irule = irule + 1;
    end
else%nrules == 0 - Unregulated Gene
    pactnr = 0;%since pact = pactnr/(1+pactdr)
    pactdr = 0;
end
pact_nr = [];
pact_dr = [];
connectlogic = connectlogic(~cellfun('isempty',connectlogic));
if ~isempty(connectlogic)
    for iconnect = 1:nrules-1
        if strcmp(connectlogic{iconnect},'|') && isempty(pact_nr)
            pact_nr = pactnr(iconnect) + pactnr(iconnect+1);
        elseif strcmp(connectlogic{iconnect},'|')
            pact_nr = pact_nr + pactnr(iconnect+1);
        end
        if strcmp(connectlogic{iconnect},'&') && isempty(pact_dr)
            pact_nr = pactnr(iconnect)*pactnr(iconnect+1);
        elseif strcmp(connectlogic{iconnect},'&')
            pact_nr = pact_nr*pactnr(iconnect+1);
        end
        if ~isempty(pact_dr)
            pact_dr = pact_dr + pactdr(iconnect+1);
        else
            pact_dr = pactdr(iconnect)+pactdr(iconnect+1);
        end
    end
else
    pact_nr = pactnr;
    pact_dr = pactdr;
end

%%
function [terms,nprots,new_baff,inh_ind,act_ind] = returnterms_activity(GeneRules)
%Determine all the proteins in the rule
orpos = strfind(GeneRules,'|');
andpos = strfind(GeneRules,'&');        
if ~isempty(orpos) && isempty(andpos)
    terms = strsplit(GeneRules,'|');
elseif ~isempty(andpos) && isempty(orpos)
    terms = strsplit(GeneRules,'&');
elseif isempty(orpos) && isempty(andpos)
    terms = {GeneRules};
end
nprots = length(terms);  
new_name = cell(nprots,1);
new_baff = zeros(nprots,1);
new_sign = zeros(nprots,1);
for iprot = 1:nprots
    tfp = strcmp(terms{iprot},prot_name);
    if any(tfp)
        new_name{iprot} = prot_name{tfp};
        new_baff(iprot) = bindaff(tfp);
        new_sign(iprot) = RS(1,tfp);                    
    end
end
inh_ind = logical(new_sign < 0);%inhibiting indices
act_ind = logical(new_sign > 0);%activating indices

end

end

