function [pact] = pactivity(trnmodel,bindaffcell,ngene)
% ngene = 17;
%Need to check effect of multiplying transcription rate to binding
%affinities associated with an '&' logic.
%October 2013
pact = zeros(ngene,1);
for igene = 1:ngene
    nrules = length(trnmodel.GeneRules{igene});
    irule = 1;      
    tfcount = 1;
    connectlogic = cell(nrules,1);
    pactnr = zeros(nrules,1);
    pactdr = zeros(nrules,1);
    if nrules > 0
        while irule <= nrules
            orpos = strfind(trnmodel.GeneRules{igene}{irule},'|');
            andpos = strfind(trnmodel.GeneRules{igene}{irule},'&');        
            [terms,xterms] = regexp(trnmodel.GeneRules{igene}{irule},'(\w+.?)(\W+?)+','tokens','split');
            if ~isemptyr(xterms)
                terms = [terms,cell(1,1)];
                terms{end}{1} = xterms{end};
            end
            nprots = length(terms);       
            if nprots == length(orpos) + length(andpos);
                if ~isempty(orpos) && isemptyr(connectlogic{irule})
                    connectlogic{irule} = trnmodel.GeneRules{igene}{irule}(orpos(end));
                    terms{end}{2} = {};
                elseif ~isempty(andpos) && isemptyr(connectlogic{irule})
                    connectlogic{irule} = trnmodel.GeneRules{igene}{irule}(andpos(end));
                    terms{end}{2} = {};
                else
                    connectlogic{irule} = {};
                end        
            end
            terms{end} = terms{end}(~cellfun('isempty',terms{end}));
        %The portion below is wrong because the logic is defined for the
        %TF following the term being checked
        %i.e if checking logic for term{iprot} then q = q for term{iprot+1}
        %To be corrected !
        %I think this has been corrected!Fingers Crossed!

            iprot = 1;
            while iprot <= nprots
                tf = strcmp(terms{iprot}{1},trnmodel.Protein);
                if iprot == 1                
                    switch trnmodel.RS(igene,tf)
                        case 1
                            pactnr(irule) = trnmodel.srate(igene,tf)*bindaffcell{igene}(tfcount);
                            %pactnr(irule) = bindaffcell{igene}(tfcount);
                            pactdr(irule) = bindaffcell{igene}(tfcount);
                            %tfcount = tfcount + 1;
                        case -1
                            pactnr(irule) = trnmodel.srate(igene,tf)*1;
                            %pactnr(irule) = 1;
                            pactdr(irule) = (bindaffcell{igene}(tfcount))^2;
                            %tfcount = tfcount + 1;
                        case 2%Dual Regulators are assumed to be positive to begin with
                            %Hence pactivity of dual regulators is same as that
                            %of activating TFs
                            pactnr(irule) = trnmodel.srate(igene,tf)*bindaffcell{igene}(tfcount);
                            %pactnr(irule) = 1;
                            pactdr(irule) = bindaffcell{igene}(tfcount);
                            %tfcount = tfcount + 1;
                    end
                    tfcount = tfcount + 1;
                    if length(terms{iprot}) > 1
                        nextlogic = terms{iprot}{2};
                    else
                        nextlogic = {};
                    end                
                elseif iprot > 1 && iprot < nprots 
                    [pactnr,pactdr] = nestedfunc(nextlogic,trnmodel.RS(igene,tf),pactnr,pactdr);
                    nextlogic = terms{iprot}{2};                
                else %iprot == nprots
                    [pactnr,pactdr] = nestedfunc(nextlogic,trnmodel.RS(igene,tf),pactnr,pactdr);
                    nextlogic = {};
                end      
                iprot = iprot + 1;
            end 
            %pact(igene) = pact
            irule = irule + 1;        
        end
    
        if ~isemptyr(connectlogic)
            nconnectlogic = length(connectlogic);
        %ncoonectlogic = nrules
        %# non empty connect logics = nrules-1
            ilogic = 1;
        %ptempnr = zeros(nconnectlogic-1,1);
        %ptemdr = zeros(nconnectlogic-1,1);
            while ilogic <= nconnectlogic - 1
            %if ilogic == 1
                %ptempnr(ilogic) = pactnr(ilogic);
                %ptempdr(ilogic) = pactdr(ilogic);
            %elseif ilogic > 1 && ilogic < nconnectlogic-1
                if connectlogic{ilogic} == '|'                   
                    pact(igene) = sum(pactnr)/(1+sum(pactdr));
                elseif connectlogic{ilogic} == '&'
                    pact(igene) = prod(pactnr)/(1+sum(pactdr));
                end
                ilogic = ilogic + 1;
            end
        else
            pact(igene) = pactnr/(1+pactdr);
        end
    else
        %basal transcription/promoter activity = 1*basal rate
        %Change to a more sensible function later
        pact(igene) = trnmodel.brate(igene)*1/(1E-15*6.023E+23);
    end
end

function [pactnr,pactdr] = nestedfunc(logic,effect,pactnr,pactdr)
    if logic == '|'
        switch effect
            case 1
                pactnr(irule) = pactnr(irule) + trnmodel.srate(igene,tf)*bindaffcell{igene}(tfcount);
                pactdr(irule) = pactdr(irule) + bindaffcell{igene}(tfcount);
            case -1
                 pactnr(irule) = pactnr(irule) + trnmodel.srate(igene,tf)*1;
                 pactdr(irule) = pactdr(irule) + (bindaffcell{igene}(tfcount))^2;
            case 2
                 pactnr(irule) = pactnr(irule) + trnmodel.srate(igene,tf)*bindaffcell{igene}(tfcount);
                 pactdr(irule) = pactdr(irule) + bindaffcell{igene}(tfcount);
        end
        tfcount = tfcount + 1;
    elseif logic == '&'
        switch effect
            case 1
                pactnr(irule) = pactnr(irule)*bindaffcell{igene}(tfcount);
                pactdr(irule) = pactdr(irule) + bindaffcell{igene}(tfcount);
            case -1
                pactnr(irule) = pactnr(irule)*1;
                pactdr(irule) = pactdr(irule) + pactdr(irule)*(bindaffcell{igene}(tfcount))^2;
            case 2
                 pactnr(irule) = pactnr(irule)*bindaffcell{igene}(tfcount);
                 pactdr(irule) = pactdr(irule) + bindaffcell{igene}(tfcount);
        end
        tfcount = tfcount + 1;
    end 
end
        

end


    

        