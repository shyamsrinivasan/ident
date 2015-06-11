%function to alter initial values of specific variables 
%(sample initial values for MC simulation)
% [y0new,pbind,savesim,status] = MCsample_initval(model,ng,var,y0old,lb,ub,nsampl)
function [y0new,pbind,savesim,status] = MCMetabolic_InitialVal(model,var,y0old,lb,ub,nsampl)
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
savesim = cell(nvar,4);
ylb = zeros(nvar,1);
yub = zeros(nvar,1);
y0new = repmat(y0old,1,nsampl);   
fprintf('Variable Name\t\tUpperBound\t\tLowerBound\n');
               
for ivar = 1:nvar
    tfm = find(strcmpi(var{ivar},model.Metabolites));
%     tfg = find(strcmpi(var{ivar},model.Gene));
%     tfr = find(strcmpi(var{ivar},model.Regulators));
    if any(tfm)%any(tfg) && ~varfl(tfg)
        varfl(tfm) = 1;
        ylb(ivar) = lb(ivar)*y0old(tfm);
        yub(ivar) = ub(ivar)*y0old(tfm);
%         if ylb(ivar) == 0 && yub(ivar)-1e-9 == 0
%             pbind(ivar) = tfg;
%         end
        fprintf('%s\t\t\t%6.8f\t\t\t%6.8f\n',var{ivar},...
                                             ylb(ivar),...
                                             yub(ivar));
        pd = makedist('Uniform','lower',ylb(ivar),'upper',yub(ivar));  
        ySample = random(pd,nsampl,1); 
%     elseif any(tfr) && ~varfl(tfr)
%         varfl(ng(1)+tfr) = 1;
%         ylb(ivar) = lb(ivar)*y0old(ng(1)+tfr);
%         yub(ivar) = ub(ivar)*y0old(ng(1)+tfr);
%         if ylb(ivar) == 0 && yub(ivar)-1e-9 == 0
%             pbind(ivar) = ng(1)+tfr;
%         end
%         fprintf('%s\t\t\t%6.8f\t\t\t%6.8f\n',var{ivar},...
%                                              ylb(ivar),...
%                                              yub(ivar));
%         pd = makedist('Uniform','lower',ylb(ivar),'upper',yub(ivar));
%         ySample = random(pd,nsampl,1); 
    else
        fprintf('%s is not available\n Aborting!!!\n',var{ivar});
        status = -1;
    end       
    for isam = 1:nsampl
        if any(tfm)
            y0new(tfm,isam) = ySample(isam);
            savesim{ivar,1} = var{ivar};
            savesim{ivar,2} = ylb(ivar);
            savesim{ivar,3} = yub(ivar);            
%         elseif any(tfr)
%             y0new(ng(1)+tfr,isam) = ySample(isam);
%             savesim{ivar,1} = var{ivar};
%             savesim{ivar,2} = ylb(ivar);
%             savesim{ivar,3} = yub(ivar);            
        end
    end
    if any(tfm)
        savesim{ivar,4} = y0new(tfm,:);
%     elseif any(tfr)
%         savesim{ivar,4} = y0new(ng(1)+tfr,:);
    end
end


for ivar = 1:nvar
    tfm = find(strcmpi(var{ivar},model.Metabolites));
%     tfg = find(strcmpi(var{ivar},model.Gene));
%     tfr = find(strcmpi(var{ivar},model.Regulators));
    if any(tfm) %&& ~varfl(tfg)
        
        y0new(tfm) = ylb+(yub-ylb)*betarnd(1.5,4.5,length(tfm),1);
        
        savesim{ivar,1} = var{ivar};
        savesim{ivar,2} = lb(ivar);
        savesim{ivar,3} = ub(ivar);
        savesim{ivar,4} = ylb;
        savesim{ivar,5} = yub;
        savesim{ivar,6} = y0new(tfm);
        
        if ylb == 0 && yub == 0
            pbind(ivar) = tfm;
        end
%     elseif any(tfr) && ~varfl(tfr)
%         ylb = lb(ivar)*y0old(ng(1)+tfr);
%         yub = ub(ivar)*y0old(ng(1)+tfr);
%         y0new(ng(1)+tfr) = ylb+(yub-ylb)*betarnd(1.5,4.5,length(tfr),1);
%         
%         savesim{ivar,1} = var{ivar};
%         savesim{ivar,2} = lb(ivar);
%         savesim{ivar,3} = ub(ivar);
%         savesim{ivar,4} = ylb;
%         savesim{ivar,5} = yub;
%         savesim{ivar,6} = y0new(ng(1)+tfr);
%         fprintf('New %s initial Value is :%6.3f\n',model.Regulators{tfr},...
%                                                    y0new(ng(1)+tfr));
%         if ylb == 0 && yub == 0
%             pbind(ivar) = ng(1)+tfr;
%         end
    else
        fprintf('%s is not available\n Aborting!!!\n',var{ivar});
        status = -1;
    end        
end
if ~any(pbind)
    pbind = [];
end
return