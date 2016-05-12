% hfig - figure handle for sub figure
% hsfig - sub figure handle on which to plot EV yield spaces
% efm - EV or EFM whose yield is to be plotted, row wise
% idx - index of reactions correspodning to columns of efm (output from
% CNA)
% -----------------------------------------------------------------------
function [allefm,yields] = printEVyieldspace(hfig,hsfig,fluxid,efm,idx,model,subid)

% convert efm to allefm (for allrxns in model.rxns)
nefm = size(efm,1);
if isfield(model,'rxns')
    nrxn = length(model.rxns);
    reacID = model.rxns;
elseif isifield(model,'reacID')
    nrxn = size(model.reacID,1);
    reacID = cellstr(model.reacID);
end
ntidx = setdiff(1:nrxn,idx);
allefm = zeros(nefm,nrxn);

for irxn = 1:nrxn
    if ismember(irxn,idx)
        allefm(:,irxn) = efm(:,idx==irxn);
    elseif ismember(irxn,ntidx)
        allefm(:,irxn) = zeros(nefm,1);
    end
end

fig_name = sprintf('Flux Envelope');
if isempty(hfig)       
    if isempty(findobj('type','figure','Name',fig_name))
        hfig = figure('Name',fig_name);     
    else
        hfig = findobj('type','figure','Name',fig_name);
    end  
end
figure(hfig);

yield = zeros(nefm,nrxn);
textstr = cell(nefm,1);
for ihs = 1:length(hsfig)
    axes(hsfig(ihs));
    set(gca,'NextPlot','add');
    for iefm = 1:nefm
        yield(iefm,:) = allefm(iefm,:)/allefm(iefm,subid);
        if any(isnan(yield(iefm,:)))
            yield(iefm,isnan(yield(iefm,:))) = 0;
        end    
        if any(isinf(yield(iefm,:)))
            yield(iefm,isinf(yield(iefm,:))) = 0;
        end
        line(yield(iefm,fluxid(1)),yield(iefm,fluxid(2)),...
             'LineStyle','none','Marker','o','MarkerEdgeColor','r',...
             'MarkerFaceColor','r','MarkerSize',6);            
        textstr{iefm} = sprintf('EV %d',iefm);        
    end
end
yields = [yield(:,fluxid(1)) yield(:,fluxid(2))];
t = text(yields(:,1)+0.02,yields(:,2)+0.02,textstr);

