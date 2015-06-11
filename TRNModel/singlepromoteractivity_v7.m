function [pact] = singlepromoteractivity_v7(Protein,srate,brate,RS,...
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
%version 7.0
% - April 01 2014 AND boolean functions changed to behave as intended
% - Minor revisions still required to eliminate the use of defparval.srate
% 
%**************************************************************************              
%Actual Code
pact = zeros(1,1);
nrules = length(GeneRules);
%irule = 1;
%tfcount = 1;
connectlogic = cell(nrules,1);
pactnr = zeros(nrules,1);
pactdr = zeros(nrules,1);
%psrate = zeros(nrules,1);
psrate_all = [];

if nrules == 1
    if ~isempty(GeneRules)
        [terms] = return_terms(GeneRules{1});
        terms{end} = terms{end}(~cellfun('isempty',terms{end}));
    %   Determine pactivitiy
        [pactnr,pactdr,psrate,pact] =...
        promoter_act(terms,srate,bindaff,defparval,RS,Protein);      
        psrate_all = [psrate_all, psrate];
        %pact = pactnr/(1+pactdr);
    end
elseif nrules > 1
    irule = 1;
    while irule <= nrules      
        [terms,nprots] = return_terms(GeneRules{irule});
        %Determine connectlogic
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

        %Determine pactivity
        [pactnr(irule),pactdr(irule),psrate] =...
        promoter_act(terms,srate,bindaff,defparval,RS,Protein);
        %psrate_all = [psrate_all, psrate];

        irule = irule + 1;
    end
else%nrules == 0 - Unregulated Gene
    if ~isempty(bindaff)
        %return bindaff;
        pact = brate(igene)*bindaff;%remove brate;
    end
end

if ~isemptyr(connectlogic)
    for krule = 1:nrules-1
        if connectlogic{krule} == '&'
            if krule == 1
                nr_all = pactnr(krule)*pactnr(krule+1);
                dr_all = pactdr(krule) + pactdr(krule+1);
            elseif krule > 1                 
                nr_all = nr_all*pactnr(krule+1);
                dr_all = dr_all + pactdr(krule+1);
            end
            nr_all = nr_all/defparval.srate;%length(psrate_all(:,krule))/sum(psrate_all(:,krule));
            %pact = (nr_all/(1+dr_all))*length(psrate_all(:,krule))/sum(psrate_all(:,krule));
                
        elseif connectlogic{krule} == '|'
            if krule == 1
                nr_all = pactnr(krule) + pactnr(krule+1);
                dr_all = pactdr(krule) + pactdr(krule+1);
            elseif krule > 1
                nr_all = nr_all + pactnr(krule+1);
                dr_all = dr_all + pactdr(krule+1);
            end             
        end
    end
    pact = nr_all/(1+dr_all); 
end

function [terms,nprots] = return_terms(GeneRules)
%Determine all the proteins in the rule
orpos = strfind(GeneRules,'|');
andpos = strfind(GeneRules,'&');        
[terms,xterms] = regexp(GeneRules,'(\w+.?)(\W+?)+','tokens','split');
%terms = terms{1}; d
%xterms = xterms{1};
if ~isemptyr(xterms)
    terms = [terms,cell(1,1)];
    terms{end}{1} = xterms{end};
end
nprots = length(terms);  
end

end




