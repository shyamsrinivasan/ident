function pvec = VmaxFBA(model,pvec,mc)

unkwid = {'pyk','mdh','pfl','g6pdh2r','ppc','me1','me2','sucdi',...
          'frd7','NADTRHD','NADH16','THD2','ATPS4r'};

vcollc = zeros(model.nt_rxn,1);      
for irxn = 1:length(unkwid)
    vid = strcmpi(model.rxns,unkwid{irxn});
    %if ismember(find(vid),model.Vind)
        [~,ck] = CKinetics(model,pvec,mc,find(vid));
    %elseif ismember(find(vid),model.Vex)
        %[~,ck] = TKinetics(model,pvec,mc,find(vid));
    %end
    if ck
        pvec.Vmax(vid) = model.Vss(vid)/ck;
    else
        pvec.Vmax(vid) = 0;
    end
    vcollc(vid) = 1;
end

%Vmax for exchange and other transport reactions
Vex = model.Vex;
for irxn = 1:length(Vex)
    if ~ismember(Vex(irxn),find(vcollc))
        [~,ck] = TKinetics(model,pvec,mc,Vex(irxn));
        if ck
            pvec.Vmax(Vex(irxn)) = model.Vss(Vex(irxn))/ck;
        elseif ~ck && model.Vss(Vex(irxn))
            pvec.Vmax(Vex(irxn)) = 1;
        elseif ~ck && ~model.Vss(Vex(irxn))
            pvec.Vmax(Vex(irxn)) = 0;
        end           
    end
end


