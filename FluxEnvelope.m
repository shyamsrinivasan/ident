function [hsubfig,prxnid,flag] = FluxEnvelope(model,varname,prodid)
if nargin < 3
    prxnid = strcmpi('Pex',model.rxns);
else
    prxnid = prodid;
end
if nargin < 2
    rxnid = setdiff(1:model.nt_rxn,model.VFup);    
else
    rxnid = [];
    for iv = 1:length(varname)
        tfp = strcmpi(varname{iv},model.rxns);
        if any(tfp)
            rxnid = union(rxnid,find(tfp));
        end
    end   
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

%Find maximum/minimum allowable growth rate
%Uptake Flux
bounds.Vuptake = model.Vuptake;
bounds.vl = zeros(model.nt_rxn,1);
% bounds.vl(bounds.vl==0) = -1;
% bounds.vl(bounds.vl==0) = -100;
bounds.vl(logical(model.rev)) = -100;%bounds.Vuptake;
bounds.vu = zeros(model.nt_rxn,1);          
%Corresponding flux bounds
bounds.vu(bounds.vu==0) = 100;%bounds.Vuptake;
[gMax,gMin,~,~,Maxflag,Minflag] = solveLP(model,'','',bounds,model.bmrxn);
if Maxflag > 0 && Minflag > 0
    fprintf('\nMaximum Allowable growth Rate = %2.3g h-1\n',-gMax);    
    fprintf('Minimum Allowable growth Rate = %2.3g h-1\n',gMin);
end
if model.gmax > -gMax
    model.gmax = -gMax;
end

%Find maximum/minimum allowable growth rate
[pMax,pMin,~,~,Maxflag,Minflag] = solveLP(model,'','',bounds,find(prxnid));
if Maxflag > 0 && Minflag > 0
    fprintf('\nMaximum allowable product flux = %2.3g h-1\n',-pMax);    
    fprintf('Minimum allowable product flux = %2.3g h-1\n',pMin);
end
    
% Max_flag = zeros(1,length(flval),model.nt_rxn);
% Min_flag = zeros(1,length(flval),model.nt_rxn);
for ifl = 1:nrxn
    %Uptake Flux
    bounds.Vuptake = model.Vuptake;
    bounds.vl = zeros(model.nt_rxn,1);
    % bounds.vl(bounds.vl==0) = -1;
%     bounds.vl(logical(model.rev)) = -100;
    bounds.vl(bounds.vl==0) = -100;   
    bounds.vu = zeros(model.nt_rxn,1);          
    %Corresponding flux bounds
    bounds.vu(bounds.vu==0) = 100;
    %Determine Max and Min for flux to be constrained with =
    [vMax,vMin,~,~,Maxflag,Minflag] = solveLP(model,'','',bounds,rxnid(ifl));
    if Maxflag > 0 && Minflag > 0
        flval(:,ifl) = linspace(vMin,-vMax,npts);
    else
        flval(:,ifl) = [];
    end
    for iv = 1:size(flval,1)
        bounds.vl(rxnid(ifl)) = flval(iv,ifl);
        bounds.vu(rxnid(ifl)) = flval(iv,ifl);        
%         prxnid = strcmpi('Pex',model.rxns);
        [fmax,fmin,vLPmax,vLPmin,Maxflag,Minflag] =...
        solveLP(model,'P','P5',bounds,find(prxnid));
%         Max_flag(1,iv,ifl) = Maxflag;
%         Min_flag(1,iv,ifl) = Minflag;
        if ~isempty(vLPmax) && Maxflag > 0
            vLPmxAll(:,iv,ifl) = vLPmax(1:model.nt_rxn);
            Maxtarget(1,iv,ifl) = fmax;            
        end
        if ~isempty(vLPmin) && Minflag > 0
            vLPmnAll(:,iv,ifl) = vLPmin(1:model.nt_rxn);
            Mintarget(1,iv,ifl) = fmin;
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
% fig_name = texlabel(['Flux Envelope \mu = ' num2str(model.gmax) 'h^{-1}']);
fig_name = sprintf('Flux Envelope mu = %g h-1',model.gmax);
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
        hsubfig(rxnid(ifl)) = subplot(nrows,2,ifl);
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