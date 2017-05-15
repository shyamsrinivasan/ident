function [model,par] = addRxn(model,par,rxn)
model.K = par.K;
model.Klb = par.Klb;
model.Kub = par.Kub;
model.KIact = par.KIact;
model.KIihb = par.KIihb;
model.Vmax = par.Vmax;
model.kcat_fwd = par.kcat_fwd;
model.kcat_bkw = par.kcat_bkw;

if isfield(rxn,'equation')
    if iscell(rxn.equation)
        eqstr = rxn.equation;
    else
        eqstr{1} = rxn.equation;
    end
end

if isfield(rxn,'par')
    if iscell(rxn.par)
        allpar = rxn.par;
    else
        allpar{1} = char(rxn.par);
    end
else
    allpar = cell(length(eqstr),1);
end

for iadd = 1:length(eqstr)
%     ipt = 0;
    imetab = length(model.mets)+1;
    irxn = length(model.rxns)+1;

    [model,reverse] =...
    GetRxnInfo(model,imetab,irxn,allpar{iadd},eqstr{iadd});    
   
    model.rev(irxn) = reverse;
    
    if isfield(rxn,'name')
        model.rxns{irxn} = rxn.name{iadd};
        model.enzs{irxn} = rxn.name{iadd};
    else
        model.rxns{irxn} = sprintf('newRxnID_%d',iadd);
        model.enzs{irxn} = sprintf('newRxnID_%d',iadd);
    end
    
    if isfield(rxn,'delSGr')
        model.delSGr(irxn) = rxn.delSGr(iadd);
    elseif isfield(rxn,'Keq')
        model.delSGr(irxn) = -0.008314*298*log(Keq(irxn));
    else
        model.delSGr(irxn) = 0;
    end
    
    if isfield(rxn,'delGlb')
        model.delGlb(irxn) = rxn.delGlb(irxn);
    else
        model.delGlb(irxn) = 0;
    end
    
    if isfield(rxn,'delGub')
        model.delGub(irxn) = rxn.delGub(irxn);
    else
        model.delGub(irxn) = 0;
    end
    
    if isfield(rxn,'Keq')
        model.Keq(irxn) = rxn.Keq(irxn);
    elseif isfield(model,'delSGr')
        model.Keq(irxn) = exp(-model.delSGr(irxn)/(0.008314*298));
    else
        model.Keq(irxn) = 0;
    end
    
    if isfield(rxn,'Vss')
        model.Vss(irxn) = rxn.Vss(irxn);
    else
        model.Vss(irxn) = 0;
    end
    
    if isfield(rxn,'kcat_fwd')
        model.kcat_fwd(irxn) = rxn.kcat_fwd(irxn);
    else
        model.kcat_fwd(irxn) = 0;
    end

    if isfield(rxn,'kcat_bkw')
        model.kcat_bkw(irxn) = rxn.kcat_bkw(irxn);
    else
        model.kcat_bkw(irxn) = 0;
    end

    if isfield(rxn,'Vmax')
        model.Vmax(irxn) = rxn.Vmax(irxn);
    else
        model.Vmax(irxn) = 0;
    end
    
    [nmet,~] = size(model.S);  
    %regulator additions do not yet work
    model.SI = [model.SI sparse(nmet,1)]; 
    model.KIact = [model.KIact sparse(nmet,1)];
    model.KIihb = [model.KIihb sparse(nmet,1)];
end

[model,par] = rearrange_par(model);
[model.nt_metab,model.nt_rxn] = size(model.S);
[nt_metab,nt_rxn] = size(model.S);    
    
model.b = zeros(nt_metab,1);
model.c = sparse(1,model.bmrxn,1,1,nt_rxn)';
model.lb = zeros(nt_rxn,1);
model.lb(model.lb==0) = -100;
model.lb(model.bmrxn) = 0;
model.ub = zeros(nt_rxn,1);
model.ub(model.ub==0) = 100; 

function [model,par] = rearrange_par(model)

%Separate External & Internal mets
exter_mind = ~cellfun('isempty',regexp(model.mets,'\w(?:\[e\])$'));
inter_mind = ~exter_mind;
bm_ind = strcmpi(model.mets,'Biomass[c]');
inter_mind(bm_ind)=0;

model.mets = [model.mets(inter_mind,1);...
              model.mets(exter_mind,1);...
              model.mets(bm_ind,1)];
model.S = [model.S(inter_mind,:);model.S(exter_mind,:);model.S(bm_ind,:)];
model.SI = [model.SI(inter_mind,:);model.SI(exter_mind,:);model.SI(bm_ind,:)];
model.K = [model.K(inter_mind,:);model.K(exter_mind,:);model.K(bm_ind,:)];
model.Kl = [model.Klb(inter_mind,:);model.Klb(exter_mind,:);model.Klb(bm_ind,:)];
model.Ku = [model.Kub(inter_mind,:);model.Kub(exter_mind,:);model.Kub(bm_ind,:)];
newKIa = [model.KIact(inter_mind,:);...
          model.KIact(exter_mind,:);...
          model.KIact(bm_ind,:)];
newKIi = [model.KIihb(inter_mind,:);...
          model.KIihb(exter_mind,:);...
          model.KIihb(bm_ind,:)];

[~,nt_rxn] = size(model.S);   
%New Indices
% model.rxns = model.rxns;
[Vind,VFex,Vex,bmrxn] = fluxIndex(model,nt_rxn,model.S);

nenz = length(model.enzs);
other_ind = setdiff(1:nenz,[Vind...                            
                            Vex...
                            VFex'...
                            bmrxn]);
mS = size(model.S,1); 
model.enzs = [model.enzs(setdiff(1:nenz,other_ind),1);...
              model.enzs(other_ind,1)];
model.rxns = [model.rxns(setdiff(1:nenz,other_ind),1);...
              model.rxns(other_ind,1)];    
          
model.S = [model.S,sparse(mS,length(other_ind))];
model.SI = [model.SI,sparse(mS,length(other_ind))];
model.K = [model.K,sparse(mS,length(other_ind))];
model.Kl = [model.Kl,sparse(mS,length(other_ind))];
model.Ku = [model.Ku,sparse(mS,length(other_ind))];
model.KIact = [newKIa,sparse(mS,length(other_ind))];
model.KIihb = [newKIi,sparse(mS,length(other_ind))];
model.Vss = [model.Vss;zeros(length(other_ind),1)];
model.delSGr = [model.delSGr;zeros(length(other_ind),1)];
model.delGlb = [model.delGlb;zeros(length(other_ind),1)];
model.delGub = [model.delGub;zeros(length(other_ind),1)];
model.Keq = [model.Keq;zeros(length(other_ind),1)];
model.kcat_fwd = [model.kcat_fwd;zeros(length(other_ind),1)];
model.kcat_bkw = [model.kcat_bkw;zeros(length(other_ind),1)];
model.rev = [model.rev;zeros(length(other_ind),1)];          

[nt_metab,nt_rxn] = size(model.S);   

%New Indices
[Vind,VFex,Vex,bmrxn] = fluxIndex(model,nt_rxn,model.S);

model.Vind = Vind;
model.Vex = Vex;
model.VFex = VFex;
model.bmrxn = bmrxn;

%Identify activated reactions
[~,allactrxn] = find(model.SI(:,1:nt_rxn)>0);
model.Vact_ind = unique(allactrxn);

%Identify Inhibited reactions
[~,allihbrxn] = find(model.SI(:,1:nt_rxn)<0);
model.Vihb_ind = unique(allihbrxn);

model.n_rxn = length(model.Vind);
model.nt_metab = nt_metab;
model.next_metab = length(find(exter_mind));
model.nint_metab = length(find(inter_mind));

par.K = model.K;
par.Klb = model.Klb;
par.Kub = model.Kub;
par.KIact = model.KIact;
par.KIihb = model.KIact;
par.Vmax = model.Vmax;
par.kcat_fwd = model.kcat_fwd;
par.kcat_bkw = model.kcat_bkw;

model = rmfield(model,{'K','Klb','Kub','KIact','KIihb','Vmax','kcat_fwd','kcat_bkw'});




