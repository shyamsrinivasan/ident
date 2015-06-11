%==========================================================================
%create a new generule when adding newer regulatory factors with a
%different logic
%==========================================================================
function [GeneRule] = changeGeneRules(GeneRule,regulator,symb)
if nargin < 3 || isempty(symb)
    symb = '|';
end
% function [GeneRule] = changeGeneRules(GeneRule,regulator)
if ~isempty(GeneRule)
    nrules = length(GeneRule);  
    nxtlogic = cell(nrules,1);
    if nrules > 0
        irule = 1;
        while irule <= nrules    
            orpos = strfind(GeneRule{irule},'|');
            andpos = strfind(GeneRule{irule},'&');
            if ~isempty(orpos)
                nprots = length(orpos) + 1;
                for iprot = 1:nprots
                    nxtlogic{irule}(iprot) = '|';
                end                
            elseif ~isempty(andpos)
                nprots = length(andpos) + 1;
                for iprot = 1:nprots
                    nxtlogic{irule}(iprot) = '&';
                end  
            elseif isempty(orpos) && isempty(andpos)
                nprots = nrules;
            end
            
%             [terms,xterms] = regexp(GeneRule{irule},'(\w+.?)(\W+?)+','tokens','split');           
%             if ~isemptyr(xterms)
%                 terms = [terms,cell(1,1)];
%                 terms{end}{1} = xterms{end};
%             end
%             nprots = length(terms);    
%             iprot = 1;
%              while iprot <= nprots          
%                 if length(terms{iprot}) > 1
%                     nxtlogic{irule}(iprot) = terms{iprot}{2};
%                 end
%                 iprot = iprot + 1;
%              end
             
             if length(nxtlogic{irule}) == nprots                 
                 irule = irule + 1;
                 continue;
             elseif length(nxtlogic{irule}) < nprots                 
                 if ~isempty(nxtlogic{irule}) &&...
                     length(nxtlogic{irule}(nxtlogic{irule} == '|')) == length(nxtlogic{irule})                
                     GeneRule{irule} = [GeneRule{irule},sprintf('|%s',regulator)]; 
                 elseif ~isempty(nxtlogic{irule}) &&...
                         length(nxtlogic{irule}(nxtlogic{irule} == '&')) == length(nxtlogic{irule})
                     GeneRule{irule} = [GeneRule{irule},sprintf('&%s',regulator)]; 
                 elseif isempty(nxtlogic{irule})
                     GeneRule{irule} = [GeneRule{irule},sprintf('%s%s',symb,regulator)];
                 else                
                    GeneRule{irule} = [GeneRule{irule},'|'];          
                    GeneRule = [GeneRule;sprintf('%s',regulator)];
                 end          
             end            
            irule = irule + 1;
        end
    end
else
    GeneRule = regulator;
end
end