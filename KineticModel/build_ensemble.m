%Sample Quantities to build ensemble of models
function [ensb,flag] = build_ensemble(nmodels,model,pmeter,MC)
% nmodels = 100;
nsamples = nmodels;

%glcpts
Mglc = strcmpi('glc[e]',model.Metabolites);
Vglc = find(model.S(Mglc,:)<0);

try
    Vind = [model.Vind;Vglc];
catch
    Vind = [model.Vind Vglc];
end
nrxn = length(Vind);
nt_rxn = model.nt_rxn;
ntmetab = size(model.S,1);
for isample = 1:nsamples
    mname = sprintf('model%d',isample);
    ensb.(mname).K = pmeter.K;
    ensb.(mname).Kind = sparse(ntmetab,nt_rxn);
    ensb.(mname).KIact = pmeter.KIact;
    ensb.(mname).KIihb = pmeter.KIihb;
    ensb.(mname).Vmax = pmeter.Vmax;
end   
subs = cell(nrxn,1);
pruds = cell(nrxn,1);
nmetab = zeros(nrxn,1);
act = cell(nrxn,1);
inhib = cell(nrxn,1);
nreg = zeros(nrxn,1);
for irxn = 1:nrxn
    subs{irxn} = find(model.S(:,Vind(irxn))<0);
    pruds{irxn} = find(model.S(:,Vind(irxn))>0);
    nmetab(irxn) = length(find(model.S(:,Vind(irxn))));
    act{irxn} = find(model.SI(:,Vind(irxn))>0);
    inhib{irxn} = find(model.SI(:,Vind(irxn))<0);
    nreg(irxn) = length(find(model.SI(:,Vind(irxn))));
end
  %Remove kinetic consideration for metabolites like water which are
  %non-limiting always
for isample = 1:nsamples
    mname = sprintf('model%d',isample);  
    jmetab = 0;
    kreg = 0;    
    for irxn = 1:nrxn   
        Kscol = zeros(ntmetab,1);
        Kpcol = zeros(ntmetab,1);
        %Substrate/Product Samples            
        met_sat = betarnd(1.5,4.5,nmetab(irxn),1); %Beta Distribution 
        nsubs = length(subs{irxn});        
        if nsubs ~= 0
            sub_ratio = met_sat(1:nsubs)./(1-met_sat(1:nsubs));
            Ksubs = MC(subs{irxn})./sub_ratio;
            Kscol(subs{irxn},1) = ensb.(mname).K(subs{irxn},Vind(irxn));
            ensb.(mname).K(Kscol==1,Vind(irxn)) =...
            Ksubs(ensb.(mname).K(subs{irxn},Vind(irxn))==1); 
            ensb.(mname).Kind(Kscol==1,Vind(irxn)) = 1;
        end
        if nmetab(irxn)-nsubs ~= 0
            prod_ratio = met_sat(nsubs+1:nmetab(irxn))./(1-met_sat(nsubs+1:nmetab(irxn)));
            Kprod = MC(pruds{irxn})./prod_ratio;
            Kpcol(pruds{irxn},1) = ensb.(mname).K(pruds{irxn},Vind(irxn));
            ensb.(mname).K(Kpcol==1,Vind(irxn)) =...
            Kprod(ensb.(mname).K(pruds{irxn},Vind(irxn))==1);
            ensb.(mname).Kind(Kpcol==1,Vind(irxn)) = 1;
        end
        %Regulator Samples             
        reg_sat = betarnd(1.5,4.5,nreg(irxn),1);
        nact = length(act{irxn});
        if nact ~= 0
            act_ratio = reg_sat(1:nact)./(1-reg_sat(1:nact));
            ensb.(mname).KIa(act{irxn},Vind(irxn)) = MC(act{irxn})./act_ratio;
%             ensb.(mname).KI(kreg+1:kreg+nact) = MC(act{irxn})./act_ratio;
        end
        if nreg(irxn)-nact ~= 0
            inhib_ratio = reg_sat(nact+1:nreg(irxn))./(1-reg_sat(nact+1:nreg(irxn))); 
            ensb.(mname).KIi(inhib{irxn},Vind(irxn)) = MC(inhib{irxn})./inhib_ratio;
%             ensb.(mname).KI(kreg+nact+1:kreg+nreg(irxn)) = MC(inhib{irxn})./inhib_ratio;
        end
        jmetab = jmetab + nmetab(irxn);
        kreg = kreg + nreg(irxn);
    end   
    [~,vflux] = ConvinienceKinetics(model,ensb.(mname),MC,model.bmrxn,Vind);
    ensb.(mname).Vmax(Vind) = model.Vss(Vind)./vflux(Vind);
    ensb.(mname).Vmax(setdiff(1:nt_rxn,Vind)) = 1;    
end

%Select models from ensemble
flag = 0;
for is = 1:nsamples
     mname = sprintf('model%d',is);  
    if any(ensb.(mname).Vmax < 0)
        flag = 1;        
    end
end