function [newVmax,feasible] = getVmax(Vmax,model,pvec,mc,Vind,fhandle)
newVmax = pvec.Vmax;
Vmax = Vmax(Vind);
rho = model.rho;

if all(pvec.check(Vind)>0)
    [~,cflx] = fhandle(model,pvec,mc,Vind);
    tempVmax = model.Vss(Vind)./(3600.*cflx).*rho;
    tempVmax(cflx==0) = 1;
    tempVmax(~isnan(Vmax)) = Vmax(~isnan(Vmax));
    newVmax(Vind) = tempVmax; 
    % zero ss flux means no enzyme for reaction in kinetic model
    newVmax(model.Vss==0) = 0;
    if ~any(newVmax(Vind)<0)
        feasible = 1;
    else
        feasible = 0;
    end
else
    feasible = 0;
end