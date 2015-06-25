%function to alter initial values of specific variables 
%(sample initial values for MC simulation)
function [y0new,pbind,savesim,status] = MCsample_initval(model,ng,var,y0old,lb,ub,nsampl)
if isempty(ng)
    mdes = 2;
else
    mdes = 1;
end
varfl = zeros(length(y0old),1);
y0new = y0old;
nvar = length(var);
status = 0;
if nvar == 1 && strcmpi(var{1},'all')
    ylb = lb*y0;
    yub = ub*y0;
    y0new = ylb+(yub-ylb).*betarnd(1.5,4.5,length(y0new),1);
    return
end
pbind = zeros(nvar,1);
savesim = cell(nvar,6);
ylb = zeros(nvar,1);
yub = zeros(nvar,1);
y0new = repmat(y0old,1,nsampl);   
fprintf('Variable Name\t\tUpperBound\t\tLowerBound\n');
               
for ivar = 1:nvar
    if mdes == 1%trnmodel
        tf1 = find(strcmpi(var{ivar},model.Gene));%tfg
        tf2 = find(strcmpi(var{ivar},model.Regulators));%tfr        
    elseif mdes == 2 %kinmodel
        tf1 = find(strcmpi(var{ivar},model.mets));
        tf2 = [];      
    end
    
    if any(tf1) && ~varfl(tf1) 
        varfl(tf1) = 1;
        ylb(ivar) = lb(ivar)*y0old(tf1);
        yub(ivar) = ub(ivar)*y0old(tf1);
        if mdes == 1
            if ylb(ivar) == 0 && yub(ivar)-1e-9 == 0
                pbind(ivar) = tf1;
            end
        end
        fprintf('%s\t\t\t%6.8f\t\t\t%6.8f\n',var{ivar},...
                                             ylb(ivar),...
                                             yub(ivar));
        pd = makedist('Uniform','lower',ylb(ivar),'upper',yub(ivar));  
        ySample = random(pd,nsampl,1); 
        
    elseif ~isempty(tf2) && any(tf2) && ~varfl(tf2)
        varfl(ng(1)+tf2) = 1;
        ylb(ivar) = lb(ivar)*y0old(ng(1)+tf2);
        yub(ivar) = ub(ivar)*y0old(ng(1)+tf2);        
        if ylb(ivar) == 0 && yub(ivar)-1e-9 == 0
            pbind(ivar) = ng(1)+tf2;
        end
        fprintf('%s\t\t\t%6.8f\t\t\t%6.8f\n',var{ivar},...
                                             ylb(ivar),...
                                             yub(ivar));
        pd = makedist('Uniform','lower',ylb(ivar),'upper',yub(ivar));
        ySample = random(pd,nsampl,1);
    else
        fprintf('%s is not available\n Aborting!!!\n',var{ivar});
        status = -1;
        return
    end    
    
    for isam = 1:nsampl
        if any(tf1)
            y0new(tf1,isam) = ySample(isam);                       
        elseif ~isempty(tf2) && any(tf2)
            y0new(ng(1)+tf2,isam) = ySample(isam);                   
        end
    end
    savesim{ivar,1} = var{ivar};
    savesim{ivar,2} = lb(ivar);
    savesim{ivar,3} = ub(ivar);
    savesim{ivar,4} = ylb(ivar);
    savesim{ivar,5} = yub(ivar); 
    if any(tf1)
        savesim{ivar,6} = y0new(tf1,:);
    elseif any(tf2)
        savesim{ivar,6} = y0new(ng(1)+tf2,:);
    end
end

if ~any(pbind)
    pbind = [];
end
return