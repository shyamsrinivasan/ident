function model = converttoSBMLformat(FBAmodel)

model.id = 'ToyN1';
model.name = 'Klamts Toy Network N1';
model.species = struct('id',FBAmodel.mets);

for jrxn = 1:FBAmodel.nt_rxn
    met_sb = FBAmodel.mets(FBAmodel.S(:,jrxn)<0);
    met_pr = FBAmodel.mets(FBAmodel.S(:,jrxn)>0);
    mets = [met_sb;met_pr];
    S_sb = -full(FBAmodel.S(FBAmodel.S(:,jrxn)<0,jrxn));
    S_pr = full(FBAmodel.S(FBAmodel.S(:,jrxn)>0,jrxn));    
    Ss = [S_sb;S_pr];
    if ~isempty(met_sb)
        sbstr = met_sb{1};
    else
        sbstr = [];
    end
    if ~isempty(met_pr)
        prstr = met_pr{1};
    else
        prstr = [];
    end
    for isb = 2:length(met_sb)
        if ~isempty(sbstr)
            if S_sb(isb)>1
                sbstr = [sbstr sprintf('+%d%s',S_sb(isb),met_sb{isb})];
            else
                sbstr = [sbstr sprintf('+%s',met_sb{isb})];
            end
        else            
            if S_sb(isb)>1
                sbstr = [sbstr sprintf('%d%s',S_sb(isb),met_sb{isb})];
            else
                sbstr = [sbstr sprintf('%s',met_sb{isb})];
            end
        end
    end
    for ipr = 2:length(met_pr)
        if ~isempty(prstr)
            if S_pr(ipr)>1
                prstr = [prstr sprintf('+%d%s',S_pr(ipr),met_pr{ipr})];
            else
                prstr = [prstr sprintf('+%s',met_pr{ipr})];
            end
        else            
            if S_pr(ipr)>1
                prstr = [prstr sprintf('%d%s',S_pr(ipr),met_pr{ipr})];
            else
                prstr = [prstr sprintf('%s',met_pr{ipr})];
            end
        end
    end
    if FBAmodel.rev(jrxn)
        rxnstr = sprintf('%s<->%s',sbstr,prstr);
    elseif ~FBAmodel.rev(jrxn)
        rxnstr = sprintf('%s->%s',sbstr,prstr);
    end
    model.reaction(jrxn) = struct('id',rxnstr,...
                                  'reactant',struct('species',met_sb),...                                                    
                                  'product',struct('species',met_pr),...                                                   
                                  'reversible',FBAmodel.rev(jrxn));
    for isb = 1:length(met_sb)
        model.reaction(jrxn).reactant(isb).stoichiometry(1) = S_sb(isb);
    end
    for ipr = 1:length(met_pr)
        model.reaction(jrxn).product(ipr).stoichiometry(1) = S_pr(ipr);
    end    
end