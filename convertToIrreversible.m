function modelIrrev = convertIrreversible(model)

%declare variables
modelIrrev.S = spalloc(size(model.S,1),0,2*nnz(model.S));
modelIrrev.rxns = [];
modelIrrev.rev = zeros(2*length(model.rxns),1);
modelIrrev.lb = zeros(2*length(model.rxns),1);
modelIrrev.ub = zeros(2*length(model.rxns),1);
modelIrrev.c = zeros(2*length(model.rxns),1);
matchRev = zeros(2*length(model.rxns),1);
 
nRxns = size(model.S,2);
irrev2rev = zeros(2*length(model.rxns),1);

%loop through each column/rxn in the S matrix building the irreversible
%model
cnt = 0;
for i = 1:nRxns
    cnt = cnt + 1;
     
    %expand the new model (same for both irrev & rev rxns
    modelIrrev.rev(cnt) = model.rev(i);
    irrev2rev(cnt) = i;

    % Check if reaction is declared as irreversible, but bounds suggest
    % reversible (i.e., having both positive and negative bounds
    if (model.ub(i) > 0 && model.lb(i) < 0) && model.rev(i) == false
         model.rev(i) = true;
        warning(cat(2,'Reaction: ',model.rxns{i},' is classified as irreversible, but bounds are positive and negative!'))

    end
    
     % Reaction entirely in the negative direction
     if (model.ub(i) <= 0 && model.lb(i) < 0)
         % Retain original bounds but reversed
         modelIrrev.ub(cnt) = -model.lb(i);
         modelIrrev.lb(cnt) = -model.ub(i);
         % Reverse sign
         modelIrrev.S(:,cnt) = -model.S(:,i);
         modelIrrev.c(cnt) = -model.c(i);
         modelIrrev.rxns{cnt} = [model.rxns{i} '_r'];
        model.rev(i) = false;
         modelIrrev.rev(cnt) = false;
     else
         % Keep positive upper bound
         modelIrrev.ub(cnt) = model.ub(i);
         %if the lb is less than zero, set the forward rxn lb to zero
         if model.lb(i) < 0
             modelIrrev.lb(cnt) = 0;
         else
             modelIrrev.lb(cnt) = model.lb(i);
         end
         modelIrrev.S(:,cnt) = model.S(:,i);
         modelIrrev.c(cnt) = model.c(i);
         modelIrrev.rxns{cnt} = model.rxns{i}; 
     end
 
    
     %if the reaction is reversible, add a new rxn to the irrev model and
     %update the names of the reactions with '_f' and '_b'
     if model.rev(i) == true
        cnt = cnt + 1;   
         modelIrrev.rxns{cnt-1} = [model.rxns{i} '_f'];
         modelIrrev.S(:,cnt) = -model.S(:,i);
         modelIrrev.rxns{cnt} = [model.rxns{i} '_b'];
         modelIrrev.rev(cnt) = true;
         modelIrrev.lb(cnt) = 0;
         modelIrrev.ub(cnt) = -model.lb(i);
         modelIrrev.c(cnt) = 0;
     end
end

% Build final structure
modelIrrev.S = modelIrrev.S(:,1:cnt);
modelIrrev.ub = columnVector(modelIrrev.ub(1:cnt));
modelIrrev.lb = columnVector(modelIrrev.lb(1:cnt));
modelIrrev.c = columnVector(modelIrrev.c(1:cnt));
modelIrrev.rev = modelIrrev.rev(1:cnt);
modelIrrev.rev = columnVector(modelIrrev.rev == 1);
modelIrrev.rxns = columnVector(modelIrrev.rxns); 
modelIrrev.mets = model.mets;
matchRev = columnVector(matchRev(1:cnt));
modelIrrev.match = matchRev;
 
modelIrrev.b = model.b;
modelIrrev.reversibleModel = false;