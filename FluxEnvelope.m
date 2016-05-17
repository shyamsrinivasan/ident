% draw production envelope between any 2 fluxes 'flux1'(x-axis - typically 
% biomass) and 'flux2'(y-axis - typically target)
% flux1 and flux2 can be a cell array of strings or just double indices 
% eg call: FluxEnvelope(FBAmodel,{'PGI','exPYR';'PFK','exPYR'},ess_rxn);

function [hfig,hsubfig,fluxid,flag] = FluxEnvelope(model,fluxid,ess_rxn)
if nargin < 3
    ess_rxn = {};
end
if nargin < 2
    error('flenlope:flux2id','Specify fluxes to draw an envelope');
else
    if size(fluxid,2)<2
        error('flenlope:flux2id','Need atleast 2 fluxes to draw an envelope');
    end
end
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
%     flux1id = cellfun(@(x)strcmpi(x,flux1id),model.rxns,'UniformOutput',false);
%     flux1id = cell2mat(cellfun(@(x)any(x),flux1id,'UniformOutput',false));
%     flux1id = find(flux1id);   
%     
%     flux2id = cellfun(@(x)strcmpi(x,flux2id),model.rxns,'UniformOutput',false);
%     flux2id = cell2mat(cellfun(@(x)any(x),flux2id,'UniformOutput',false));
%     flux2id = find(flux2id);
%     
%     if length(flux2id)<length(flux1id)
%         flux2id = repmat(flux2id,length(flux1id),1);
%     end
%     if length(flux1id)<length(flux2id)
%         flux1id = repmat(flux1id,length(flux2id),1);
%     end
elseif ~iscell(fluxid)
    flux1id = fluxid(:,1);
    flux2id = fluxid(:,2);
end

% number of different envelopes (x(flux1)-y(flux2) combinations)
nrxn = size(fluxid,1);  
nt_rxn = size(model.S,2);
flag = 1;
npts = 100;

vLPmxAll = zeros(nt_rxn,npts,nrxn);
vLPmnAll = zeros(nt_rxn,npts,nrxn);
Maxtarget = zeros(1,npts,nrxn);
Mintarget = zeros(1,npts,nrxn);
flval = zeros(npts,nrxn);

for irxn = 1:nrxn
    
    % fix flux bounds as in FBA for uptake fluxes
    [model,bounds] = changebounds(model,ess_rxn);
    
    if bounds.vl(flux1id(irxn))==bounds.vu(flux1id(irxn))
        if bounds.vu(flux1id(irxn)) > 0
            bounds.vl(flux1id(irxn))= 0;
        elseif bounds.vu(flux1id(irxn)) == 0
            bounds.vl(flux1id(irxn))= -1;
        end
    end
    % Find maximum/minimum allowable flux 1 (x-axis for envelope)
    [LP1max,LP1min] = solveLP(model,bounds,ess_rxn,flux1id(irxn));
    if LP1max.flag > 0 && LP1min.flag > 0
        fprintf('\nMaximum Allowable %s flux = %2.3g h-1\n',...
            model.rxns{flux1id(irxn)},-LP1max.obj);    
        fprintf('Minimum Allowable %s flux = %2.3g h-1\n',...
            model.rxns{flux1id(irxn)},LP1min.obj);
    end
    
    if bounds.vl(flux2id(irxn))==bounds.vu(flux2id(irxn))
        if bounds.vu(flux2id(irxn)) > 0
            bounds.vl(flux2id(irxn))= 0;
        elseif bounds.vu(flux2id(irxn)) == 0
            bounds.vl(flux2id(irxn))= -1;
        end
    end
    % Find maximum/minimum allowable flux 2 (y-axis for envelope)
    [LP2max,LP2min] = solveLP(model,bounds,ess_rxn,flux2id(irxn));
    if LP2max.flag > 0 && LP2min.flag > 0
        fprintf('\nMaximum allowable %s flux = %2.3g h-1\n',...
            model.rxns{flux2id(irxn)},-LP2max.obj);    
        fprintf('Minimum allowable %s flux = %2.3g h-1\n',...
            model.rxns{flux2id(irxn)},LP2min.obj);
    end
    
    % Determine Max and Min for flux to be constrained with =
    if LP1max.flag > 0 && LP1min.flag > 0
        flval(:,irxn) = linspace(LP1min.obj,-LP1max.obj,npts);
    else
        flval(:,irxn) = [];
    end
    
    % fix flux1 at each of flval(ivalue) and max/min flux2
    for iv = 1:npts
        bounds.vl(flux1id(irxn)) = flval(iv,irxn);
        bounds.vu(flux1id(irxn)) = flval(iv,irxn);        

        [LPmax,LPmin] = solveLP(model,bounds,ess_rxn,flux2id(irxn));
        if ~isempty(LPmax.v) && LPmax.flag > 0
            vLPmxAll(:,iv,irxn) = LPmax.v(1:nt_rxn);
            Maxtarget(1,iv,irxn) = -LPmax.obj;            
        end
        if ~isempty(LPmin.v) && LPmin.flag > 0
            vLPmnAll(:,iv,irxn) = LPmin.v(1:nt_rxn);
            Mintarget(1,iv,irxn) = LPmin.obj;
        end
    end
end

if iscell(fluxid)
    fluxid = [flux1id flux2id];
end
    
%Infeasible Problem
if ~any(Maxtarget)
    hsubfig = [];    
    flag = -1;
    return
end

%Plot Envelope
if rem(nrxn,2) == 0
    nplots = nrxn;
else
    nplots = nrxn+1;
end
nrows = nplots/2;
if nrxn > 1
    ncol = 2;
else
    ncol=1;
end
% fig_name = texlabel(['Flux Envelope \mu = ' num2str(model.gmax) 'h^{-1}']);
fig_name = sprintf('Flux Envelope');
if isempty(findobj('type','figure','Name',fig_name))
    hfig = figure('Name',fig_name); 
    figure(hfig);
else
    hfig = findobj('type','figure','Name',fig_name);
    figure(hfig);
end
% hsubfig = zeros(nt_rxn,1);
hfig = figure;
hsubfig = zeros(nrxn,1);
ifl = 1; % 1,2,3...,nplots

while ifl <= nrxn
    if hsubfig(ifl) ~= 0
        hca = findobj(hsubfig(ifl),'type','axes');
        set(hfig,'CurrentAxes',hca);  
    else % subplot is unassigned
        hsubfig(ifl) = subplot(nrows,ncol,ifl);
        hca = gca;
        % Make sure more plots can be added at end of loop
        set(hca,'NextPlot','add');     
        set(hsubfig(ifl),'NextPlot','add');
    end    
    ylabel = sprintf('Flux %s \n mmole/mmole uptake',model.rxns{flux2id(ifl)});
    hline = plot([flval(:,ifl)' fliplr(flval(:,ifl)')],...
             [Maxtarget(1,:,ifl) fliplr(Mintarget(1,:,ifl))]);
    xlabel = sprintf('Flux %s \n mmole/mmole uptake',model.rxns{flux1id(ifl)});
    set(hline,'LineWidt',2,...
              'Color',[0 0 0]);
    set(get(gca,'YLabel'),'String',ylabel);
    set(get(gca,'XLabel'),'String',xlabel);
    axis tight;   
    ifl = ifl + 1;
end

%Annotation TextBox
% set(hfig,'CurrentAxes',hsubfig(flux1id(1)));
% h = annotation('textbox',[.3 .9 .1 .1],...
%                'String',{['$\mu$ = ' num2str(model.gmax) '$h^{-1}$'],...
%                          ['Vuptake = ' num2str(model.Vuptake) ' mmole/gDCW.h']},...
%                 'Color',[0 0 0],...
%                 'BackgroundColor',[1 1 1]);

return