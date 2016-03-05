function vLPmin = run_pFBA(model,ess_rxn,Vup_struct)
% convert to irreversible model
[modelIrrev,~,~,irrev2rev]=convertIrreversible(model);

%adding sum of all fluxes to Irrev S
%add meatabolite to S
modelIrrev.S(end+1,:) = ones(size(modelIrrev.S(1,:)));
modelIrrev.b(end+1) = 0;
modelIrrev.mets{end+1} = 'fluxMeasure';

nm = size(modelIrrev.mets,1);

%add reaction to S matrix
modelIrrev.S(:,end+1) = sparse(find(strcmpi('fluxMeasure',modelIrrev.mets)),1,-1,nm,1);
modelIrrev.rxns{end+1} = 'netFlux';
modelIrrev.vl(end+1) = 0;
modelIrrev.vu(end+1) = Inf;
modelIrrev.rev(end+1) = 0;

prxnid = find(strcmpi('netFlux',modelIrrev.rxns));
bounds.vl = modelIrrev.vl;
bounds.vu = modelIrrev.vu;
[~,vLPmin] = solveLP(modelIrrev,bounds,ess_rxn,prxnid,Vup_struct);

if vLPmin.flag > 0
    modelIrrev.vl(prxnid) = vLPmin.obj;
    modelIrrev.vu(prxnid) = vLPmin.obj;
end

bounds.vl = modelIrrev.vl;
bounds.vu = modelIrrev.vu;
[~,vLPmin] = solveLP(modelIrrev,bounds,ess_rxn,prxnid,Vup_struct);
vLPmin.v(vLPmin.v<1e-7)=0;

Vss = zeros(length(model.rxns),1);
if vLPmin.flag > 0
    ir = 1;
    %reassign irreversible fluxes to original model
    while ir<=(length(modelIrrev.rxns)-1)
        if ir<(length(modelIrrev.rxns)-1)
            if irrev2rev(ir)==irrev2rev(ir+1);
%                 irrev2rev(ir)
                Vss(irrev2rev(ir))=vLPmin.v(ir)-vLPmin.v(ir+1);
                ir = ir+2;
            else
%                 irrev2rev(ir)
                Vss(irrev2rev(ir)) = vLPmin.v(ir);
                ir = ir+1;
            end
        elseif ir==(length(modelIrrev.rxns)-1)
            Vss(irrev2rev(ir)) = vLPmin.v(ir);
            ir = ir+1;
        end
    end
else
    error('run_pFBA:pFBAfeas','The pFBA problem was infeasible');
end

Vss(abs(Vss)<1e-6) = 0;
vLPmin.v = Vss;
% model.Vss = Vss;


    