function viol = delGaconsistent(model,variable)
Vind = model.Vind;
MC = variable.MC;
MClow = model.MClow;
MChigh = model.MChigh;
viol = zeros(model.nt_rxn,1);
% for irxn = 1:length(Vind)
%     sub = model.S(:,Vind(irxn))<0;
%     prud = model.S(:,Vind(irxn))>0;
%     if any(sub) && any(prud)
%         viol(Vind(irxn)) = prod(MC(sub)) -...
%                            prod(MC(prud))/model.Keq(Vind(irxn));
%     end
%     if viol(Vind(irxn))*model.Vss(Vind(irxn)) > 0
%         viol(Vind(irxn)) = 1;
%     else
%         viol(Vind(irxn)) = -1;
%     end
% end
for irxn = 1:length(Vind)
    fprintf('Reaction #%d\n',irxn);
    flag = 1;
    while flag
        sub = model.S(:,Vind(irxn))<0;
        prud = model.S(:,Vind(irxn))>0;
        if any(sub) && any(prud)
            viol(Vind(irxn)) = prod(MC(sub)) -...
                               prod(MC(prud))/model.Keq(Vind(irxn));
                           
            if viol(Vind(irxn))*model.Vss(Vind(irxn)) >= 0
                flag = 0;
            else
                pd = makedist('Uniform');
%                 MC(sub) = 1+(10-1).*random(pd,length(find(sub)),1);
%                 MC(prud) = 1+(10-1).*random(pd,length(find(prud)),1);
                MC(sub) = MClow(sub) + (MChigh(sub) - MClow(sub)).*...
                            random(pd,length(find(sub)),1);
                MC(prud) = MClow(prud) + (MChigh(prud) - MClow(prud)).*...
                            random(pd,length(find(prud)),1);
            end
        end
        
    end
end


MC;
return