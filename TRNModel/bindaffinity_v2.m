%getbindingaffinity.m
%igene - index passed into the function 
%RS,GeneRules
function [bindaff] = bindaffinity_v2(Protein,Coefficient,RS,GeneRules,igene,protconc)
bindaff = zeros(length(find(RS(1,:))),1);
nregterm = length(GeneRules);
 
iregterm = 1;
newterms = {};    
while iregterm <= nregterm
    [terms,xterms] = regexp(GeneRules{iregterm},'(\w+.?)(\W+?)+','tokens','split');
    
    if ~isemptyr(xterms)
        terms = [terms,cell(1,1)];
        terms{end}{1} = xterms{end};
    end
    newterms = [newterms,terms];
    iregterm = iregterm + 1;
end

nterms = length(newterms);
kterm = 1;
while kterm <= nterms
    tf = strcmp(newterms{kterm}{1},Protein);
    if any(tf)            
        bindaff(kterm) = protconc(tf)/Coefficient(igene,tf);
%         if igene > 1
%             bindaff(kterm) = protconc(tf)/Coefficient(find(tf)+(igene-1)*size(tf,1),1);
%         else%if igene == 1
%             bindaff(kterm) = protconc(tf)/Coefficient(tf,1);
%         end        
        
        kterm = kterm + 1;
    end
    %kterm = kterm + 1;
end 
end