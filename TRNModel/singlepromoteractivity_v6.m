function [pact] = singlepromoteractivity_v6(Protein,srate,brate,RS,...
                  GeneRules,bindaff,igene,defparval)
%**************************************************************************
%Calculates promoter activity for genes regulated by multipe TFs over a
%single promoter
%******Input(s)
%Protein -   List of proteins from trnmodel.Protein
%srate -     Vector of stimulated transcription rates for a single gene
%            corresponding to regulating TFs
%brate -     Basal/Leakage transcriptional rate for constitutively
%            expressed/repressed genes
%RS -        Vector of regulatory effects for a single gene 
%GeneRules - Genetic regulatory rules for a single gene in cell array
%format
%bindaff -   Binding affinity of regulators regulating a gene obtained from
%            bindaffinity.m
%igene -     Gene index 
%defparval - Structure of user specified or otherwise default parameter 
%            values 
%******Output(s)
%pact - Promoter activity as a function TF concentration, binding 
%       coefficients and transcription rates  
% 
%November 2013
%version 1.0
% - Cannot incorporate custom hill coefficients - fixed at 2
% - Incorrect pactivity function inclusive of srate (TF basis)
%version 2.0
% - Incorporates srate on a gene basis as opposed to a TF basis 
%version 3.0
% - Correct pactivity function (from v2.0) for gene basis srate
%version 4.0
% - pactivity for a TF based srate function
% - December 20 2013 No corrections made for AND boolean functions yet
%version 6.0
% - March 24 2014 Correction made for AND boolean functions
% - Model failed to converge even using a stiff solver in SUNDIALS
%**************************************************************************
pact = zeros(1,1);
nrules = length(GeneRules);
irule = 1;
tfcount = 1;
connectlogic = cell(nrules,1);
pactnr = zeros(nrules,1);
pactdr = zeros(nrules,1);
pactsr = zeros(nrules,1);

if nrules > 0    
    nxtlogic = cell(nrules,1);    
    while irule <= nrules
        orpos = strfind(GeneRules{irule},'|');
        andpos = strfind(GeneRules{irule},'&');        
        [terms,xterms] = regexp(GeneRules{irule},'(\w+.?)(\W+?)+','tokens','split');
        %terms = terms{1}; d
        %xterms = xterms{1};
        if ~isemptyr(xterms)
            terms = [terms,cell(1,1)];
            terms{end}{1} = xterms{end};
        end
        nprots = length(terms);       
        if nprots == length(orpos) + length(andpos);
            if ~isempty(orpos) && isemptyr(connectlogic{irule})
                connectlogic{irule} = GeneRules{irule}(orpos(end));
                terms{end}{2} = {};
            elseif ~isempty(andpos) && isemptyr(connectlogic{irule})
                connectlogic{irule} = GeneRules{irule}(andpos(end));
                terms{end}{2} = {};
            else
                connectlogic{irule} = {};
            end        
        end
        terms{end} = terms{end}(~cellfun('isempty',terms{end}));    
         
        iprot = 1;
        nr_term = zeros(nprots,1);
        dr_term = zeros(nprots,1);
        sr_term = zeros(nprots,1);
        effect = zeros(nprots,1);
        
        while iprot <= nprots          
            if length(terms{iprot}) > 1
                nxtlogic{irule}(iprot) = terms{iprot}{2};%Try this in Debug mode!
            end
            tf = strcmp(terms{iprot}{1},Protein);
            switch RS(1,tf)
                case 1
                    nr_term(iprot) = bindaff(tfcount);
                    dr_term(iprot) = bindaff(tfcount);  
                case -1
                    nr_term(iprot) = 1;
                    dr_term(iprot) = (bindaff(tfcount))^defparval.rephill;
                case 2
                    nr_term(iprot) = bindaff(tfcount);
                    dr_term(iprot) = bindaff(tfcount);  
            end
            effect(iprot) = RS(1,tf);
            if srate(1,tf) > 0
                sr_term(iprot) = srate(1,tf);
            else
                sr_term(iprot) = defparval.srate;
            end
            tfcount = tfcount + 1;
            iprot = iprot + 1;
        end
        
        %Check to see if all logic are the same in nxtlogic
        logic_flag = 0;
        if length(nxtlogic{irule}(nxtlogic{irule} == '|')) == length(nxtlogic{irule}) 
            rule_logic = '|';
            logic_flag = 1;
        elseif length(nxtlogic{irule}(nxtlogic{irule} == '&')) == length(nxtlogic{irule})
            rule_logic = '&';
            logic_flag = 1;
        else
            rule_logic = '|';%Contigency in case
            fprintf('All Logical Conditions within a Rule should be similar \n');
            fprintf('Split Clusters with similar conditions into separate rules \n');
            %Determine separate clusters?
            %Then proceed with the usual method?      
        end
        
        if ~isemptyr(connectlogic)
            if connectlogic{irule} == '&'
                if logic_flag && rule_logic == '|'
                    pactnr(irule) = prod(nr_term(effect == -1))*(sum(nr_term(effect == 1))+sum(nr_term(effect == 2)));
                elseif logic_flag && rule_logic == '&'
                    pactnr(irule) = prod(nr_term);
                end                    
                pactsr(irule) = sum(sr_term)/length(sr_term);
            elseif connectlogic{irule} == '|'
                if logic_flag && rule_logic == '|'
                    pactnr(irule) = prod(nr_term(effect == -1))*sum(sr_term(effect == -1))/length(sr_term(effect == -1)) +...
                        sum(nr_term(effect == 1).*sr_term(effect == 1)) + ...
                        sum(nr_term(effect == 2).*sr_term(effect == 2));
                elseif logic_flag && rule_logic == '&'
                    pactnr(irule) = prod(nr_term)*sum(sr_term)/length(sr_term);
                end 
            end
        else %No connectlogic or connectlogic isempty = true
            if logic_flag && rule_logic == '|'
                pactnr(irule) = prod(nr_term(effect==-1))*sum(sr_term(effect == -1))/length(sr_term(effect == -1)) +...
                        sum(nr_term(effect == 1).*sr_term(effect == 1)) +...
                        sum(nr_term(effect == 2).*sr_term(effect == 2));
            elseif logic_flag && rule_logic == '&'
                pactnr(irule) = prod(nr_term)*sum(sr_term(sr_term > 0))/length(sr_term);
            end
        end
        pactdr(irule) = sum(dr_term);        
        irule = irule + 1;
    end
    
    if ~isemptyr(connectlogic)
        for iconnect = 1:length(connectlogic)
            if connectlogic{iconnect} == '&'
                pact = prod(pactnr)*sum(pactsr)/(length(pactsr)*(1 + sum(pactdr)));
            elseif connectlogic{iconnect} == '|'
                pact = sum(pactnr)/(1 + sum(pactdr));
            end
        end
    else
        pact = pactnr/(1 + pactdr);
    end
end
         
end