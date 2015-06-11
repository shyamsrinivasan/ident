%function [pact] = singlepromoteractivity_v2(Protein,srate,brate,RS,...
%                                         GeneRules,bindaff,igene,defparval)
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
%**************************************************************************
function [pact] = singlepromoteractivity_v4(Protein,srate,brate,RS,...
                  GeneRules,bindaff,igene,defparval)

% if nargin < 8              
%     defhillcoeff = 2;              
% end
              
pact = zeros(1,1);
nrules = length(GeneRules);
irule = 1;
tfcount = 1;
connectlogic = cell(nrules,1);
pactnr = zeros(nrules,1);
pactdr = zeros(nrules,1);
psrate = zeros(nrules,1);%New statement for psrate - stores srates for each
                         %rule. srates could be overwritten 
%pactsrate = zeros(1,1);
%defsrate = 1.66e-6;
if nrules > 0
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
        psrate_prot = zeros(nprots,1);%New statement for psrate1 - stores all 
                                  %srates individually w/o overwritting
        %old_tf = [];
        iprot = 1;
        old_tf = zeros(length(Protein),nprots);
        while iprot <= nprots
            tf = strcmp(terms{iprot}{1},Protein);
            %=========First Protein/TF in the Gene Regulatory Rule=========
            if iprot == 1             
                switch RS(1,tf)
                    case 1
                        pactnr(irule) = srate(1,tf)*bindaff(tfcount);                        
%                         pactnr(irule) = bindaff(tfcount);
                        pactdr(irule) = bindaff(tfcount);  
%                         psrate(irule) = srate(1,tf);
                    case -1
%                         pactnr(irule) = srate(igene,tf)*1;
                        pactnr(irule) = 1;
                        pactdr(irule) = (bindaff(tfcount))^defparval.rephill;
%                         psrate(irule) = srate(1,tf);
                    case 2 %pactivity of dual regulators is same as +ve TFs
                        pactnr(irule) = srate(1,tf)*bindaff(tfcount);
%                         pactnr(irule) = bindaff(tfcount);
                        pactdr(irule) = bindaff(tfcount);
%                         psrate(irule) = srate(1,tf);
                end            
                %Store this tf.
                %old_tf = [old_tf,tf];
                old_tf(:,iprot) = tf;
                %prev_tf = tf;
                tfcount = tfcount + 1;
                if length(terms{iprot}) > 1
                    nextlogic = terms{iprot}{2};
                else
                    nextlogic = {};
                end                
            elseif iprot > 1 && iprot < nprots 
                %==============All but first or last regulator=============
                [pactnr,pactdr,psrate] = nestedfunc(nextlogic,RS,iprot,pactnr,pactdr,psrate);
                %old_tf = [old_tf,tf];
                old_tf(:,iprot) = tf;
                %prev_tf = tf;
                nextlogic = terms{iprot}{2};                
            else %iprot == nprots
                %========================Last Regulator====================
                [pactnr,pactdr,psrate] = nestedfunc(nextlogic,RS,iprot,pactnr,pactdr,psrate);
                %old_tf = [old_tf,tf];
                old_tf(:,iprot) = tf;
                %prev_tf = tf;
                nextlogic = {};
            end 
            
            %Store this srate
            psrate_prot(iprot) = srate(1,tf);
            iprot = iprot + 1;
            
        end 
        if pactnr(irule) == 1
            pactnr(irule) = defparval.srate*pactnr(irule);
        end
        irule = irule + 1;        
    end
    
    if ~isemptyr(connectlogic)
        nconnectlogic = length(connectlogic);
        %ncoonectlogic = nrules
        %# non empty connect logics = nrules-1
        ilogic = 1;       
        while ilogic <= nconnectlogic - 1           
            if connectlogic{ilogic} == '|' 
                pact = sum(pactnr)/(1+sum(pactdr));   %v1
                %pact = srate*sum(pactnr)/(1+sum(pactdr)); %v2                
            elseif connectlogic{ilogic} == '&'
                %pact = prod(pactnr)/(1+sum(pactdr));  %v1
                %Wrong notation
                %Need to change entire pattern for rules invovling AND
                pact = defparval.srate*prod(pactnr)/(1+sum(pactdr)); %v2
            end
            ilogic = ilogic + 1;
        end
    else
        pact = pactnr/(1+pactdr);  %v1
        %pact = srate*pactnr/(1+pactdr); %v2
    end
else
    %basal transcription/promoter activity = 1*basal rate
    %Change to a more sensible function later
    pact = brate(igene)*1/(1E-15*6.023E+23);
end

%==================Function to calculate promoter activity=================

function [pactnr,pactdr,psrate] = nestedfunc(logic,RS,kprot,pactnr,pactdr,psrate)
%     effect = RS(1,tf);
    prev_effect = RS(1,logical(old_tf(:,kprot-1)));
    if logic == '|'
        switch RS(1,tf)
            case 1
                if prev_effect == 1
                    pactnr(irule) = pactnr(irule) + srate(1,tf)*bindaff(tfcount);
                elseif prev_effect == -1 %previous_effect = -1
                    if pactnr(irule) == 1
                        pactnr(irule) = pactnr(irule)*srate(1,tf)*bindaff(tfcount);
                    else
                        pactnr(irule) = pactnr(irule) + srate(1,tf)*bindaff(tfcount);
                    end
                end                
                pactdr(irule) = pactdr(irule) + bindaff(tfcount);
            case -1
                if prev_effect == 1
                    pactnr(irule) = pactnr(irule)*1;                    
                elseif prev_effect == -1 %previous_effect = -1
                    if pactnr(irule) == 1
                        pactnr(irule) = pactnr(irule)*1;
                    else
                        pactnr(irule) = pactnr(irule)*1;
                    end
                end
                pactdr(irule) = pactdr(irule) + (bindaff(tfcount))^defparval.rephill;               
            case 2
                if prev_effect == 1
                    pactnr(irule) = pactnr(irule) + srate(1,tf)*bindaff(tfcount);
                elseif prev_effect == -1 %previous_effect = -1
                    if pactnr(irule) == 1
                        pactnr(irule) = pactnr(irule)*srate(1,tf)*bindaff(tfcount);
                    else
                        pactnr(irule) = pactnr(irule) + srate(1,tf)*bindaff(tfcount);
                    end
                end                
                pactdr(irule) = pactdr(irule) + bindaff(tfcount);                
        end       
    elseif logic == '&'
        switch RS(1,tf)
            case 1
                if prev_effect == 1
                    pactnr(irule) = pactnr(irule)*bindaff(tfcount);
                elseif prev_effect == -1
                    if pactnr(irule) == 1
                        pactnr(irule) = pactnr(irule)*bindaff(tfcount);
                    else
                        pactnr(irule) = pactnr(irule)*bindaff(tfcount);
                    end
                end
                pactdr(irule) = pactdr(irule) + bindaff(tfcount);
            case -1
                if prev_effect == 1
                    pactnr(irule) = pactnr(irule)*1;
                elseif prev_effect == -1
                    if pactnr(irule) == 1
                        pactnr(irule) = pactnr(irule)*1;
                    else
                        pactnr(irule) = pactnr(irule)*1;
                    end
                end
                pactdr(irule) = pactdr(irule) + pactdr(irule)*(bindaff(tfcount))^defparval.rephill;
            case 2
                if prev_effect == 1
                    pactnr(irule) = pactnr(irule)*bindaff(tfcount);
                elseif prev_effect == -1
                    if pactnr(irule) == 1
                        pactnr(irule) = pactnr(irule)*bindaff(tfcount);
                    else
                        pactnr(irule) = pactnr(irule)*bindaff(tfcount);
                    end
                end
                pactdr(irule) = pactdr(irule) + bindaff(tfcount);                 
        end        
    end 
    tfcount = tfcount + 1;
end
end
        

