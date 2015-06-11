clc
% v0 = vLP;
load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Kinetic Model\N2m.mat');
bounds.vl = zeros(FBAmodel.nt_rxn,1);
% bounds.vl(bounds.vl==0) = -1;
bounds.vl(5) = -1;
bounds.vu = zeros(FBAmodel.nt_rxn,1);
%Uptake Flux
bounds.Vuptake = 1;
%Corresponding flux bounds
bounds.vu(bounds.vu==0) = 1;

Vind = FBAmodel.Vind;
flval = linspace(0,bounds.Vuptake);

vLPmxAll = zeros(FBAmodel.nt_rxn,length(flval),FBAmodel.nt_rxn+1);
vLPmnAll = zeros(FBAmodel.nt_rxn,length(flval),FBAmodel.nt_rxn+1);
Maxtarget = zeros(1,length(flval),FBAmodel.nt_rxn+1);
Mintarget = zeros(1,length(flval),FBAmodel.nt_rxn+1);
for ifl = 2:FBAmodel.nt_rxn
    bounds.vl = zeros(FBAmodel.nt_rxn,1);
    % bounds.vl(bounds.vl==0) = -1;
    bounds.vl(5) = -1;
    bounds.vu = zeros(FBAmodel.nt_rxn,1);
    %Uptake Flux
    bounds.Vuptake = 1;
    %Corresponding flux bounds
    bounds.vu(bounds.vu==0) = 1;
    for iv = 1:length(flval)
        bounds.vl(ifl) = flval(iv);
        bounds.vu(ifl) = flval(iv);        
        prxnid = 8;
        [vLPmax,vLPmin,fmax,fmin] = solveLP(FBAmodel,'P','P5',bounds,prxnid);
        if ~isempty(vLPmax)
            vLPmxAll(:,iv,ifl) = vLPmax;
            Maxtarget(1,iv,ifl) = fmax;
        end
        if ~isempty(vLPmin)
            vLPmnAll(:,iv,ifl) = vLPmin;
            Mintarget(1,iv,ifl) = fmin;
        end
        vLPmax = [];
        vLPmin = [];
        fmax = [];
        fmin = [];
    end
end
%reversible v5 bounds
bounds.vl = zeros(FBAmodel.nt_rxn,1);
% bounds.vl(bounds.vl==0) = -1;
bounds.vl(5) = -1;
bounds.vu = zeros(FBAmodel.nt_rxn,1);
%Uptake Flux
bounds.Vuptake = 1;
%Corresponding flux bounds
bounds.vu(bounds.vu==0) = 1;
ifl = 5;
for iv = 1:length(flval)
    bounds.vl(ifl) = -flval(iv);
    bounds.vu(ifl) = -flval(iv);        
    prxnid = 8;
    [vLPmax,vLPmin,fmax,fmin] = solveLP(FBAmodel,'P','P5',bounds,prxnid);
    if ~isempty(vLPmax)
        vLPmxAll(:,iv,FBAmodel.nt_rxn+1) = vLPmax;
        Maxtarget(1,iv,FBAmodel.nt_rxn+1) = fmax;
    end
    if ~isempty(vLPmin)
        vLPmnAll(:,iv,FBAmodel.nt_rxn+1) = vLPmin;
        Mintarget(1,iv,FBAmodel.nt_rxn+1) = fmin;
    end
    vLPmax = [];
    vLPmin = [];
    fmax = [];
    fmin = [];
end

%Plot Envelope
if rem(FBAmodel.nt_rxn+1,2) == 0
    nplots = FBAmodel.nt_rxn+1;
else
    nplots = FBAmodel.nt_rxn+1+1;
end
nrows = nplots/2;
figure
for ifl = 2:(FBAmodel.nt_rxn+1)
    subplot(nrows,2,ifl-1);
    ylabel = sprintf('Flux %s',FBAmodel.Enzyme{prxnid});
    if ifl == 11
        hline = plot([-flval fliplr(-flval)],...
                 [-Maxtarget(1,:,ifl) fliplr(Mintarget(1,:,ifl))],...
                 'LineWidth',2,...
                 'Color',[0 .5 0]);
        xlabel = sprintf('Flux %s',FBAmodel.Enzyme{5});
    else
        hline = plot([flval fliplr(flval)],...
                 [-Maxtarget(1,:,ifl) fliplr(Mintarget(1,:,ifl))],...
                 'LineWidth',2,...
                 'Color',[0 .5 0]);
        xlabel = sprintf('Flux %s',FBAmodel.Enzyme{ifl});
    end
    
    set(get(gca,'YLabel'),'String',ylabel);
    set(get(gca,'XLabel'),'String',xlabel);
    axis tight;   
end
