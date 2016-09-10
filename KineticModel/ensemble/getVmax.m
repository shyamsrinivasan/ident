function [newVmax,feasible] = getVmax(Vmax,model,pvec,mc,Vind,fhandle)
newVmax = pvec.Vmax;

if all(pvec.check(Vind)>0)
    [~,cflx] = fhandle(model,pvec,mc,Vind);
    tempVmax = model.Vss(Vind)./(3600.*cflx);
    tempVmax(cflx==0) = 1;
    tempVmax(~isnan(Vmax(Vind))) = Vmax(~isnan(Vmax(Vind)));
    newVmax(Vind) = tempVmax;    
    feasible = 1;
else
    feasible = 0;
end