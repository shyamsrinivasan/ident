function Nimc = changeInitialCondition(model,initval,VMCpos,VMCneg,allMets)
%change intial value by the given percentage in VMC
if nargin<5
    allMets = 0;
end
if nargin<4
    VMCneg = struct([]);
end
if nargin < 3
    fprintf('Nothing to change\n');
    Nimc = initval;
    return
end

nintmet = model.nint_metab;
Nimc = initval;
if ~isempty(VMCneg)
    [mc_lb,mc_ub] = iconcentration(model,VMCneg);
    zro_met = find(initval==0);
    if any(mc_lb(zro_met))
        idx = zro_met(logical(mc_lb(zro_met)));
        fprintf('%s is already zero. \n %s cannot be below zero',model.mets(idx));
        error('initcond:ivalChk','Normalized Initial conditions cannot be negative');
    end
    if isempty(mc_ub)
        Nimc = Nimc-Nimc.*mc_lb/100;
        Nimc(initval==0) = 0;
    elseif any(mc_lb>0) && any(mc_ub>0)
        mcid = find(mc_lb);
        for id = 1:length(mcid)
            rnd_dist =...
            random(makedist('Uniform','lower',mc_lb(mcid(id))/100,...
                                      'upper',mc_ub(mcid(id))/100),...
                                      length(mc_lb),1);
            Nimc(mcid(id)) = Nimc(mcid(id)).*rnd_dist(mcid(id));
        end
    end
end

if ~isempty(VMCpos)
    [mc_lb,mc_ub] = iconcentration(model,VMCpos);
    if isempty(mc_ub)
        Nimc = Nimc+Nimc.*mc_lb/100;
    elseif any(mc_lb>0) && any(mc_ub>0)
        mcid = find(mc_lb);
        for id = 1:length(mcid)
            rnd_dist =...
            random(makedist('Uniform','lower',mc_lb(mcid(id))/100,...
                                      'upper',mc_ub(mcid(id))/100),...
                                      1,1);
            Nimc(mcid(id)) = Nimc(mcid(id)).*rnd_dist;
        end        
    end    
end
%Nimc = Nimc + Nimc*random percentage(0-1)*Scaling factor
    %use vector direction determination from ACHR sampler for
    %concentrations
if allMets
    pd = makedist('Uniform','lower',-10,'upper',10);
    [mc_lb,mc_ub] = iconcentration(model,[]);
    if isempty(mc_ub)
        nt_sign = sign(random(pd,length(mc_lb),1));    
        rnd_dist = random(makedist('Uniform'),length(mc_lb),1);
        Nimc(1:nintmet) = Nimc(1:nintmet) +...
        nt_sign(1:nintmet).*Nimc(1:nintmet).*rnd_dist(1:nintmet)*0.5;
    elseif ~isempty(mc_lb) && ~isempty(mc_ub)
        nt_sign = sign(random(pd,length(mc_lb),1));
        if any(mc_lb>0) && any(mc_ub>0)
            rnd_dist = random(makedist('Uniform','lower',mc_lb/100,'upper',mc_ub/100),length(mc_lb),1);
            Nimc(1:nintmet) = Nimc(1:nintmet).*rnd_dist(1:nintmet);
        elseif ~any(mc_lb>0) && ~any(mc_ub>0)
            rnd_dist = random(makedist('Uniform'),length(mc_lb),1);
            Nimc(1:nintmet) = Nimc(1:nintmet) +...
            nt_sign(1:nintmet).*Nimc(1:nintmet).*rnd_dist(1:nintmet);
        end
    end        
end