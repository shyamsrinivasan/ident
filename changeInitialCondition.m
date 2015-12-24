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

Nimc = initval;
if ~isempty(VMCneg)
    mc = iconcentration(model,VMCneg);
    zro_met = find(initval==0);
    if any(mc(zro_met))
        idx = zro_met(logical(mc(zro_met)));
        fprintf('%s is already zero. \n %s cannot be below zero',model.mets(idx));
        error('initcond:ivalChk','Normalized Initial conditions cannot be negative');
    end
%     if any(mc(initval==0)
%         
    Nimc = Nimc-Nimc.*mc/100;
    Nimc(initval==0) = 0;
end

if ~isempty(VMCpos)
    mc = iconcentration(model,VMCpos);
    Nimc = Nimc+Nimc.*mc/100;
elseif allMets    
    %normed vector spaces(scalar)
    %note: vector*scalar = vector with different length but same direction
    mc = iconcentration(model,[]);
    rnd_dist = random(makedist('Uniform'),length(mc),1);
    Nimc = Nimc + Nimc.*rnd_dist;
    %Nimc = Nimc + Nimc*random percentage(0-1)*Scaling factor
    %use vector direction determination from ACHR sampler for
    %concentrations
end