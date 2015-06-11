% function [pact] = singlepromoteractivity_v7(Protein,srate,brate,RS,...
%                   GeneRules,bindaff,igene,defparval)
%Pseudo code

%while irule <= nrules
%if rule > 1
%   use new code
%else
%   use old code written as a function
%end
%end

%Change Pseudo Code to
%if nrules > 1
%   while irule <= nrules
%else
%Statements w/o loop for irules

%Actual Code
pact = zeros(1,1);
nrules = length(GeneRules);
irule = 1;
tfcount = 1;
connectlogic = cell(nrules,1);
pactnr = zeros(nrules,1);
pactdr = zeros(nrules,1);
psrate = zeros(nrules,1);
psrate_all = [];

if nrules > 1
  irule = 1;
  while irule <= nrules      
      %[terms,nprots] = return_terms();
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
      [pactnr(irule),pactdr(irule),psrate] = pact_copyv5(terms,srate,bindaff,defparval,RS,Protein);
      psrate_all = [psrate_all, psrate];

      irule = irule + 1;
  end
else
    %[terms,nprots] = return_terms();
    terms{end} = terms{end}(~cellfun('isempty',terms{end}));
    
%   Determine pactivitiy
    [pactnr,pactdr,psrate,pact] = pact_copyv5(terms,srate,bindaff,defparval,RS,Protein);      
    psrate_all = [psrate_all, psrate];
    pact = pactnr/(1+pactdr);

end



%==============Old Setting - Serves No Purpose==============
% - Old Code
%%% - Modified version used above

% % % while irule <= nrules
    %Determine all the proteins in the rule
% % %     orpos = strfind(GeneRules{irule},'|');
% % %     andpos = strfind(GeneRules{irule},'&');        
% % %     [terms,xterms] = regexp(GeneRules{irule},'(\w+.?)(\W+?)+','tokens','split');
% % %     %terms = terms{1}; d
% % %     %xterms = xterms{1};
% % %     if ~isemptyr(xterms)
% % %         terms = [terms,cell(1,1)];
% % %         terms{end}{1} = xterms{end};
% % %     end
% % %     nprots = length(terms);    
    %Determine if there is more than 1 rule
% % %     if nrules > 1
% % %         rule_logic = char();
        %Determine connect logic
% % %         if nprots == length(orpos) + length(andpos);
% % %             if ~isempty(orpos) && isemptyr(connectlogic{irule})
% % %                 connectlogic{irule} = GeneRules{irule}(orpos(end));
% % %                 terms{end}{2} = {};
% % %             elseif ~isempty(andpos) && isemptyr(connectlogic{irule})
% % %                 connectlogic{irule} = GeneRules{irule}(andpos(end));
% % %                 terms{end}{2} = {};
% % %             else
% % %                 connectlogic{irule} = {};
% % %             end        
% % %         end
% % %         terms{end} = terms{end}(~cellfun('isempty',terms{end}));
        
%         while iprot <= nprots
%             tf = strcmp(terms{iprot}{1},Protein);
%             if iprot == 1
%                 if length(terms{iprot}) > 1
%                     rule_logic(iprot) = terms{iprot}{2};            
%                 end
%             elseif ~isempty(rule_logic) && iprot > 1
%                 if length(terms{iprot}) > 1
%                     if strcmp(rule_logic(iprot-1),terms{iprot}{2})
%                         l_flag = 1;
%                     end
%                 end
%             end            
%             iprot = iprot + 1;
%         end 
        
%         if length(nxtlogic{irule}(nxtlogic{irule} == '|')) == length(nxtlogic{irule}) 
%             rule_logic = '|';
%             logic_flag = 1;
%         elseif length(nxtlogic{irule}(nxtlogic{irule} == '&')) == length(nxtlogic{irule})
%             rule_logic = '&';
%             logic_flag = 1;
%         else
%             rule_logic = '|';%Contigency in case
%             fprintf('All Logical Conditions within a Rule should be similar \n');
%             fprintf('Split Clusters with similar conditions into separate rules \n');
%             %Determine separate clusters?
%             %Then proceed with the usual method?      
%         end
        
        %Determine Pactivity using Method from v6
        %Should a method similar to v5 be used instead?Yes => Discard v6
% % %         [pactnr(irule),pactdr(irule),psrate] = pact_copyv5(terms,srate,bindaff,defparval,RS,Protein);
% % %         psrate_all = [psrate_all, psrate];
   
% % %     else
% % %         terms{end} = terms{end}(~cellfun('isempty',terms{end}));
        
%         while iprot <= nprots
%             tf = strcmp(terms{iprot}{1},Protein);
%             if iprot == 1
%                 if length(terms{iprot}) > 1
%                     rule_logic(iprot) = terms{iprot}{2};            
%                 end
%             elseif ~isempty(rule_logic) && iprot > 1
%                 if length(terms{iprot}) > 1
%                     if strcmp(rule_logic(iprot-1),terms{iprot}{2})
%                         l_flag = 1;
%                     end
%                 end
%             end            
%             iprot = iprot + 1;
%         end
        
        %Check to see if all logic are the same in nxtlogic
%         if length(nxtlogic{irule}(nxtlogic{irule} == '|')) == length(nxtlogic{irule}) 
%             rule_logic = '|';            
%         elseif length(nxtlogic{irule}(nxtlogic{irule} == '&')) == length(nxtlogic{irule})
%             rule_logic = '&';            
%         else
%             rule_logic = '|';%Contigency in case
%             fprintf('All Logical Conditions within a Rule should be similar \n');
%             fprintf('Split Clusters with similar conditions into separate rules \n');
%             %Determine separate clusters?
%             %Then proceed with the usual method?      
%         end
        
        %Determine pactivity differently
% % %         [pactnr,pactdr,psrate,pact] = pact_copyv5(terms,srate,bindaff,defparval,RS,Protein);      
% % %         psrate_all = [psrate_all, psrate];
% % %         pact = pactnr/(1+pactdr);
        
  
% % %     end
% % %     irule = irule + 1;
% % % end
%==============Old Setting - End========================

%         if ~isempty(connectlogic)
%             if connectlogic == '&'
%                 %follow method in v5 no matter what the individual rule
%                 %logics are
%                 [pactnr,pactdr,psrate] = pact_copyv5(terms,srate,bindaff,defparval,RS,Protein);
%                 %new nr = old nr * new bindaff*average of all srates??Yes                            
%             elseif connectlogic == '|'
%                 %Rule by Rule pact determination similar to case w/ one
%                 %rule
%                 %for irule <= nrules
%                 %pact(irule) = pact_copyv5();
%                 [pactnr,pactdr,psrate] = pact_copyv5(terms,srate,bindaff,defparval,RS,Protein);
%                 %total nr = nr rule 1+...+nr rule n
%             end
%         end    

    
if ~isempty(connectlogic)
    for krule = 1:nrules-1
        if connectlogic{krule} == '&'
            if krule == 1
                nr_all = pactnr(krule)*pactnr(krule+1);
                dr_all = pactdr(krule) + pacdr(krule+1);
            elseif krule > 1                 
                nr_all = nr_all*pactnr(krule+1);
                dr_all = dr_all + pactdr(krule+1);
            end
            nr_all = nr_all*length(psrate_all(:,krule))/sum(psrate_all(:,krule));
            %pact = (nr_all/(1+dr_all))*length(psrate_all(:,krule))/sum(psrate_all(:,krule));
                
        elseif connectlogic{krule} == '|'
            if krule == 1
                nr_all = pactnr(krule) + pactnr(krule+1);
                dr_all = pactdr(krule) + pacdr(krule+1);
            elseif krule > 1
                nr_all = nr_all + pactnr(krule+1);
                dr_all = dr_all + pactdr(krule+1);
            end
            
                    
        end
    end
    pact = nr_all/(1+dr_all); 
end

% function [terms,nprots] = return_terms(GeneRules)
Determine all the proteins in the rule
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
% end

% end




