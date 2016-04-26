
function [hsubfig,prxnid,flag] = FluxEnvelope(model,flux1id,flux2id,ess_rxn,varname)
if nargin < 5
    rxnid = 1:model.nt_rxn;    
else   
    rxnid = cellfun(@(x)strcmpi(x,varname),model.rxns,'UniformOutput',false);
    rxnid = cell2mat(cellfun(@(x)any(x),rxnid,'UniformOutput',false));
    rxnid = find(rxnid);    
end
if nargin < 4
    ess_rxn = {};
end
if nargin < 3
    error('flenlope:flux2id','Need atleast 2 fluxes 2 draw an envelope');
end

if iscell(flux1id)
    flux1id = find(strcmpi(model.rxns,flux1id{1}));
end
if iscell(flux2id)
    flux2id = find(strcmpi(model.rxns,flux2id{1}));
end

nrxn = length(rxnid);
    
flag = 1;
% v0 = vLP;
% load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Kinetic Model\N2m.mat');

npts = 100;
% rxnid = setdiff(1:model.nt_rxn,[model.bmrxn model.Vupind]);

vLPmxAll = zeros(model.nt_rxn,npts,nrxn);
vLPmnAll = zeros(model.nt_rxn,npts,nrxn);
Maxtarget = zeros(1,npts,nrxn);
Mintarget = zeros(1,npts,nrxn);
flval = zeros(npts,nrxn);

% fix flux bounds as in FBA for uptake fluxes
[model,bounds] = changebounds(model,ess_rxn);

% Find maximum/minimum allowable flux 1 (x-axis for envelope)
[LP1max,LP1min] = solveLP(model,bounds,ess_rxn,flux1id);
if LP1max.flag > 0 && LP1min.flag > 0
    fprintf('\nMaximum Allowable %s flux = %2.3g h-1\n',model.rxns{flux1id},-LP1max.obj);    
    fprintf('Minimum Allowable %s flux = %2.3g h-1\n',model.rxns{flux1id},LP1min.obj);
end
% if model.gmax > -gMax
%     model.gmax = -gMax;
% end

% Find maximum/minimum allowable flux 2 (y-axis for envelope)
[LP2max,LP2min] = solveLP(model,bounds,ess_rxn,flux2id);
if LP2max.flag > 0 && LP2min.flag > 0
    fprintf('\nMaximum allowable %s flux = %2.3g h-1\n',model.rxns{flux2id},-LP2max.obj);    
    fprintf('Minimum allowable %s flux = %2.3g h-1\n',model.rxns{flux2id},LP2min.obj);
end
    
for ifl = 1:nrxn
    %Uptake Flux
%     bounds.Vuptake = model.Vuptake;
%     bounds.vl = zeros(model.nt_rxn,1);
%     % bounds.vl(bounds.vl==0) = -1;
% %     bounds.vl(logical(model.rev)) = -100;
%     bounds.vl(bounds.vl==0) = -100;   
%     bounds.vu = zeros(model.nt_rxn,1);          
%     %Corresponding flux bounds
%     bounds.vu(bounds.vu==0) = 100;
    % Determine Max and Min for flux to be constrained with =
    [LPmax,LPmin] = solveLP(model,bounds,ess_rxn,rxnid(ifl));
    if LPmax.flag > 0 && LPmin.flag > 0
        flval(:,ifl) = linspace(LPmin.obj,-LPmax.obj,npts);
    else
        flval(:,ifl) = [];
    end
    for iv = 1:size(flval,1)
        bounds.vl(rxnid(ifl)) = flval(iv,ifl);
        bounds.vu(rxnid(ifl)) = flval(iv,ifl);        

        [LPmax,LPmin] = solveLP(model,bounds,ess_rxn,flux2id);
        if ~isempty(LPmax.v) && LPmax.flag > 0
            vLPmxAll(:,iv,ifl) = LPmax.v(1:model.nt_rxn);
            Maxtarget(1,iv,ifl) = -LPmax.obj;            
        end
        if ~isempty(LPmin.v) && LPmin.flag > 0
            vLPmnAll(:,iv,ifl) = LPmin.v(1:model.nt_rxn);
            Mintarget(1,iv,ifl) = LPmin.obj;
        end
        vLPmax = [];
        vLPmin = [];
        fmax = [];
        fmin = [];
    end
end

%Infeasible Problem
if ~any(Maxtarget)
    hsubfig = [];
    prxnid = [];
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
if length(varname)>1
    ncol = 2;
else
    ncol=1;
end
% fig_name = texlabel(['Flux Envelope \mu = ' num2str(model.gmax) 'h^{-1}']);
fig_name = sprintf('Flux Envelope %g',model.gmax);
if isempty(findobj('type','figure','Name',fig_name))
    hfig = figure('Name',fig_name); 
    figure(hfig);
end
hsubfig = zeros(model.nt_rxn,1);
ifl = 1;%1,2,3...,nplots
%trxid = rxid(ifl);%1,2,3,...,nt_rxn
while ifl <= nrxn
% for ifl = 2:(nrxn)
    if hsubfig(rxnid(ifl)) ~= 0
%     if hsubfig(ifl) ~= 0
        hca = findobj(hsubfig(rxnid(ifl)),'type','axes');
%         hca = findobj(hsubfig(ifl),'type','axes');
        set(hfig,'CurrentAxes',hca);  
    else%subplot is unassigned
        hsubfig(rxnid(ifl)) = subplot(nrows,ncol,ifl);
%         hsubfig(ifl) = subplot(nrows,2,ifl-1);   
        hca = gca;
        %Make sure more plots can be added at end of loop
        set(hca,'NextPlot','add');     
        set(hsubfig(rxnid(ifl)),'NextPlot','add');
    end    
    ylabel = sprintf('Flux %s \n mmole/mmole uptake',model.rxns{prxnid});
%     if ifl == 11
%         hline = plot([-flval fliplr(-flval)]./model.Vuptake,...
%                  [-Maxtarget(1,:,ifl) fliplr(Mintarget(1,:,ifl))]./model.Vuptake);
%         xlabel = sprintf('Flux %s',model.rxns{5});
%     else
        hline = plot([flval(:,ifl)' fliplr(flval(:,ifl)')],...
                 [-Maxtarget(1,:,ifl) fliplr(Mintarget(1,:,ifl))]);
        xlabel = sprintf('Flux %s \n mmole/mmole uptake',model.rxns{rxnid(ifl)});
%         xlabel = sprintf('Flux %s',model.rxns{ifl});
%     end
    set(hline,'LineWidt',2,...
              'Color',[0 0 0]);
    set(get(gca,'YLabel'),'String',ylabel);
    set(get(gca,'XLabel'),'String',xlabel);
    axis tight;   
    ifl = ifl + 1;
end

%Annotation TextBox
set(hfig,'CurrentAxes',hsubfig(rxnid(1)));
h = annotation('textbox',[.3 .9 .1 .1],...
               'String',{['$\mu$ = ' num2str(model.gmax) '$h^{-1}$'],...
                         ['Vuptake = ' num2str(model.Vuptake) ' mmole/gDCW.h']},...
                'Color',[0 0 0],...
                'BackgroundColor',[1 1 1]);

return