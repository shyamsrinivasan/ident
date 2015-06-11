% function [pact] = singlepromoteractivity(Protein,srate,brate,RS,...
%                                          GeneRules,bindaff,igene)
%**************************************************************************
%version 1.0
% - Cannot incorporate custom hill coefficients - fixed at 2
% - Incorrect pactivity function inclusive of srate (TF basis)
%**************************************************************************
function [pact] = singlepromoteractivity(Protein,srate,brate,RS,GeneRules,bindaff,igene)
pact = zeros(1,1);
nrules = length(GeneRules);
irule = 1;
tfcount = 1;
connectlogic = cell(nrules,1);
pactnr = zeros(nrules,1);
pactdr = zeros(nrules,1);
if nrules > 0
    while irule <= nrules
        orpos = strfind(GeneRules{irule},'|');
        andpos = strfind(GeneRules{irule},'&');        
        [terms,xterms] = regexp(GeneRules{irule},'(\w+.?)(\W+?)+','tokens','split');
        %terms = terms{1};
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
        while iprot <= nprots
            tf = strcmp(terms{iprot}{1},Protein);
            if iprot == 1                
                switch RS(1,tf)
                    case 1
                        pactnr(irule) = srate(igene,tf)*bindaff(tfcount);
                        pactdr(irule) = bindaff(tfcount);                  
                    case -1
                        pactnr(irule) = srate(igene,tf)*1;
                        pactdr(irule) = (bindaff(tfcount))^2;                        
                    case 2 %pactivity of dual regulators is same as +ve TFs
                        pactnr(irule) = srate(igene,tf)*bindaff(tfcount);
                        pactdr(irule) = bindaff(tfcount);                        
                end
                tfcount = tfcount + 1;
                if length(terms{iprot}) > 1
                    nextlogic = terms{iprot}{2};
                else
                    nextlogic = {};
                end                
            elseif iprot > 1 && iprot < nprots 
                [pactnr,pactdr] = nestedfunc(nextlogic,RS(1,tf),pactnr,pactdr);
                nextlogic = terms{iprot}{2};                
            else %iprot == nprots
                [pactnr,pactdr] = nestedfunc(nextlogic,RS(1,tf),pactnr,pactdr);
                nextlogic = {};
            end      
                iprot = iprot + 1;
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
            pact = sum(pactnr)/(1+sum(pactdr));
        elseif connectlogic{ilogic} == '&'
            pact = prod(pactnr)/(1+sum(pactdr));
        end
        ilogic = ilogic + 1;
    end
    else
    pact = pactnr/(1+pactdr);
    end
else
    %basal transcription/promoter activity = 1*basal rate
    %Change to a more sensible function later
    pact = brate(igene)*1/(1E-15*6.023E+23);
end

function [pactnr,pactdr] = nestedfunc(logic,effect,pactnr,pactdr)
    if logic == '|'
        switch effect
            case 1
                pactnr(irule) = pactnr(irule) + srate(igene,tf)*bindaff(tfcount);
                pactdr(irule) = pactdr(irule) + bindaff(tfcount);
            case -1
                pactnr(irule) = pactnr(irule)*1; 
%                 pactnr(irule) = pactnr(irule) + srate(igene,tf)*1;
                pactdr(irule) = pactdr(irule) + (bindaff(tfcount))^2;
            case 2
                pactnr(irule) = pactnr(irule) + srate(igene,tf)*bindaff(tfcount);
                pactdr(irule) = pactdr(irule) + bindaff(tfcount);
        end
        tfcount = tfcount + 1;
    elseif logic == '&'
        switch effect
            case 1
                pactnr(irule) = pactnr(irule)*bindaff(tfcount);
                pactdr(irule) = pactdr(irule) + bindaff(tfcount);
            case -1
                pactnr(irule) = pactnr(irule)*1;
                pactdr(irule) = pactdr(irule) + pactdr(irule)*(bindaff(tfcount))^2;
            case 2
                 pactnr(irule) = pactnr(irule)*bindaff(tfcount);
                 pactdr(irule) = pactdr(irule) + bindaff(tfcount);
        end
        tfcount = tfcount + 1;
    end 
end
end
        

