function plotPointsonFluxEnvelope(hfig,hsubfig,fluxid,flux)

if nargin < 4
    error('flptlope:data','Specify data to be plotted');
end
if nargin < 3
    error('flptlope:flux2id','Specify fluxes to draw an envelope');
else
    if size(fluxid,2)<2
        error('flptlope:flux2id','Need atleast 2 fluxes to draw an envelope');
    end
end
if iscell(fluxid)
    flux1id = fluxid(:,1);    
    flux1id = cellfun(@(x)strcmpi(x,flux1id),model.rxns,'UniformOutput',false);
    flux1id = cell2mat(cellfun(@(x)any(x),flux1id,'UniformOutput',false));
    flux1id = find(flux1id);
    
    flux2id = fluxid(:,2);
    flux2id = cellfun(@(x)strcmpi(x,flux2id),model.rxns,'UniformOutput',false);
    flux2id = cell2mat(cellfun(@(x)any(x),flux2id,'UniformOutput',false));
    flux2id = find(flux2id);
    
    if length(flux2id)<length(flux1id)
        flux2id = repmat(flux2id,length(flux1id),1);
    end
    if length(flux1id)<length(flux2id)
        flux1id = repmat(flux1id,length(flux2id),1);
    end
elseif ~iscell(fluxid)
    flux1id = fluxid(:,1);
    flux2id = fluxid(:,2);
end

if size(fluxid,1) ~= length(hsubfig)
    error('flptlope:szmatch',...
    'Size mismatch between number of available figures and number of flux IDs');
end

nrxn = size(fluxid,1);
for ifl = 1:nrxn
    set(hfig,'CurrentAxes',hsubfig(ifl));
    set(hsubfig(ifl),'NextPlot','add');
    
    for iss = 1:size(flux,2)
        cid = chooseColors(1,{'blue'});
        line(flux(flux1id(ifl),iss),flux(flux2id(ifl),iss),'LineStyle','none',...
                             'Marker','.','MarkerEdgeColor',cid{1},...
                             'MarkerFaceColor',cid{1},'MarkerSize',30);
        ptid = ['SS' num2str(iss)];
        ht = text(flux(flux1id(ifl),iss),flux(flux2id(ifl),iss),ptid,...
             'HorizontalAlignment','right',...
             'VerticalAlignment','bottom','Interpreter','latex');
         set(ht,'FontSize',22,'FontName','Arial',...
             'FontWeight','bold');
    end
    axesP.FontName  = 'Arial';
    axesP.FontSize = 22;
    axesP.LineWidth = 1.5;
    axesP.TickLength = [0.01 0.01];
    axesP.XColor = [.1 .1 .1];
    axesP.YColor = [.1 .1 .1];
    set(gca,axesP);
end


