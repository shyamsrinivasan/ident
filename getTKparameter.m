function [] = getTKparameter(model,pvec,mc)
Vex = model.Vex;

for irxn = 1:length(Vex)
    if model.Vss(Vex(irxn))~=0
        %kinetics - substrate and product
        sbid = logical(model.S(:,Vex(irxn))<0);
        prid = logical(model.S(:,Vex(irxn))>0);
        if any(sbid) && any(prid)
            if ~strcmpi(model.rxns{Vex(irxn)},'h2ot')
                sbid(h2o) = 0;
            end
            if ~strcmpi(model.rxns{Vex(irxn)},'PIt2r')
                sbid([pic pie hc he]) = 0;
                prid([pic pie hc he]) = 0;
            else
                sbid([hc he]) = 0;
                prid([hc he]) = 0;
            end
            [x,fval,exitflag,output,lambda] =...
            fmincon(@vflux,x0,[],[],[],[],lb,ub,nonlncon,options);
        end
    end
end

function flux = vflux(x)
kcatfwd = x(1);
kcatbkw = x(2);

if model.rev(Vex(irx))
flux = (kcatfwd*mc(sbid)-kcatbkw*mc(prid))./...
        (1+mc(sbid)+mc(prid));
elseif ~model.rev(Vex(irxn))
    flux = (kcatfwd*mc(sbid))./...
        (1+mc(sbid));
end
end
end

    