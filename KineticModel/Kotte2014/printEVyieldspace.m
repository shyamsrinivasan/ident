% hfig - figure handle for sub figure
% hsfig - sub figure handle on which to plot EV yield spaces
% efm - EV or EFM whose yield is to be plotted, row wise
% idx - index of reactions correspodning to columns of efm (output from
% CNA)
% efm - EFM from CNA
% idx - index of reactions in efm from CNA
% -----------------------------------------------------------------------
function [allefm,yields] = printEVyieldspace(hfig,hsfig,fluxid,efm,idx,model,subid)

if iscell(fluxid)
    nfig = size(fluxid,1);
    flux1id = fluxid(:,1);    
    flux2id = fluxid(:,2);
    for ifig = 1:nfig
        flux1id{ifig} = find(strcmpi(model.rxns,flux1id{ifig}));
        flux2id{ifig} = find(strcmpi(model.rxns,flux2id{ifig}));
    end
    flux1id = cell2mat(flux1id);
    flux2id = cell2mat(flux2id);
elseif ~iscell(fluxid)
    flux1id = fluxid(:,1);
    flux2id = fluxid(:,2);
end

% convert efm to allefm (for allrxns in model.rxns)
nefm = size(efm,1);
if isfield(model,'rxns')
    nrxn = length(model.rxns);
%     reacID = model.rxns;
elseif isifield(model,'reacID')
    nrxn = size(model.reacID,1);
%     reacID = cellstr(model.reacID);
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

if size(fluxid,1) ~= length(hsfig)
    error('EFMyield:szmatch',...
    'Size mismatch between number of available figures and number of flux IDs');
end
% figure axes properties set at the end
axesP.FontName  = 'Arial';
axesP.FontSize = 22;
axesP.LineWidth = 1.5;
axesP.TickLength = [0.01 0.01];
axesP.XColor = [.1 .1 .1];
axesP.YColor = [.1 .1 .1];

nfig = size(flux1id,1);
yield = zeros(nefm,nrxn);
textstr = cell(nefm,1);
for ifl = 1:nfig
    set(hfig,'CurrentAxes',hsfig(ifl));
    set(hsfig(ifl),'NextPlot','add');
    for iefm = 1:nefm
        yield(iefm,:) = allefm(iefm,:)/allefm(iefm,subid);
        if any(isnan(yield(iefm,:)))
            yield(iefm,isnan(yield(iefm,:))) = 0;
        end    
        if any(isinf(yield(iefm,:)))
            yield(iefm,isinf(yield(iefm,:))) = 0;
        end
        line(yield(iefm,flux1id(ifl)),yield(iefm,flux2id(ifl)),...
             'LineStyle','none','Marker','.','MarkerEdgeColor','r',...
             'MarkerFaceColor','r','MarkerSize',25); 
        ptid = ['EV' num2str(iefm)];
        ht = text(yield(iefm,flux1id(ifl)),yield(iefm,flux2id(ifl)),ptid,...
             'HorizontalAlignment','right',...
             'VerticalAlignment','baseline','Interpreter','latex');
        set(ht,'FontSize',22,'FontName','Arial');
%         textstr{iefm} = sprintf('EV %d',iefm);
         
        set(get(gca,'YLabel'),'FontName','Arial');   
        set(get(gca,'YLabel'),'FontSize',22); 
        
        set(get(gca,'XLabel'),'FontName','Arial');   
        set(get(gca,'XLabel'),'FontSize',22);
        axis tight;   
        set(gca,axesP);
    end
end


% fig_name = sprintf('Flux Envelope');
% if isempty(hfig)       
%     if isempty(findobj('type','figure','Name',fig_name))
%         hfig = figure('Name',fig_name);     
%     else
%         hfig = findobj('type','figure','Name',fig_name);
%     end  
% end
% figure(hfig);
% 
yields = [yield(:,flux1id)';yield(:,flux2id)'];
% t = text(yields(:,1)+0.02,yields(:,2)+0.02,textstr);

