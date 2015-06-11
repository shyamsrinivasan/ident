function mconc = sampleMetabolites(model,Vind)
if nargin < 2
    Vind = 1:model.n_rxn;
    n_rxn = model.n_rxn;
else
    n_rxn = length(Vind);
end
%Sample metabolites under thermodynamic constraints
% n_rxn = model.n_rxn;
nmetab = model.nint_metab;
SGr = model.delSGr/(0.008314*298);
qmA = zeros(n_rxn,1);
mconc = zeros(nmetab,1);
accpt_mind = zeros(nmetab,1);
for irxn = 1:n_rxn
    subsind = model.S(:,Vind(irxn)) < 0;
    prodind = model.S(:,Vind(irxn)) > 0;    
    s_subs = model.S(subsind,Vind(irxn));
    s_prod = model.S(prodind,Vind(irxn));
    %Sample Metabolite Concentrations
    metab_ind = [find(subsind);find(prodind)];    
    if ~all(accpt_mind(metab_ind))
        new_mind = metab_ind(~accpt_mind(metab_ind));        
        [mconc,qmA(Vind(irxn)),accpt_mind] =...
        recursive_call(mconc,model,accpt_mind,new_mind,subsind,prodind,...
        s_subs,s_prod,SGr(Vind(irxn)),model.Vss(Vind(irxn)),model.Keq(Vind(irxn)));        
    else
        %All Metabolites for the rxn are already sampled
    end
end
return

function [mconc,qmA,accpt_mind] = recursive_call(mconc,model,accpt_mind,m_ind,subsind,prodind,s_subs,s_prod,SGr,Vss,Keq)    
    mc_sampl = gen_newsample(m_ind,model);
    mconc(m_ind) = mc_sampl;
    %Mass Action Ratio
    qmA = prod(mconc(subsind).^s_subs)*prod(mconc(prodind).^s_prod);
    delGr = SGr + log(qmA);
    gamma = prod(mconc(subsind))-prod(mconc(prodind))/Keq;
    %Thermodynamic Constraint
    if Vss*gamma > 0
%     if Vss*delGr < 0
        %accept sample
        mconc(m_ind) = mc_sampl;
        accpt_mind(m_ind) = 1;
    else
        %reject sample and resample
        [mconc,qmA,accpt_mind] =... 
        recursive_call(mconc,model,accpt_mind,m_ind,subsind,prodind,s_subs,s_prod,SGr,Vss,Keq);          
    end
return

function c_sampl = gen_newsample(m_indx,model)
  nmetab_rxn = length(m_indx);
  c_sampl = model.MClow(m_indx) +...
       (model.MChigh(m_indx) - model.MClow(m_indx)).*betarnd(1.5,4.5,nmetab_rxn,1);
return



    
    
    
 

    
    
