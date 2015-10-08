nrxns = size(model.S,2);
npoints = 10;

warmup = sparse(nrxns,npoints);

i = 1;
while i <=npoints/2
    
    %create random objective function
    model.c = rnad(nrxns,1)-0.5;
    
    for maxMin = [1 -1]
        %set objective function
        if i<=nrxns
            model.c = sparse(i,1,1,nrxns,1);
        end
        model.osense=maxMin;
        
        %determine max or min
        sol = solverCobraLP(model)
        x = sol.full;
        
        %if optimal solution is found
        x(x>model.ub) = model.ub(x>model.ub);
        x(x<model.lb) = model.lb(x<model.lb);
        
        %strore points
        
        