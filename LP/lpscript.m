%developing the phase plane for the toy network
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
% % prodind = strcmpi('P5',FBAmodel.Metabolites);
% % prxnid = FBAmodel.Vexind(FBAmodel.S(prodind,FBAmodel.Vexind)~=0);
flval = linspace(0,bounds.Vuptake);
vLPmxAll = cell(FBAmodel.nt_rxn,1);
vLPmnAll = cell(FBAmodel.nt_rxn,1);
Maxtarget = cell(FBAmodel.nt_rxn,1);
Mintarget = cell(FBAmodel.nt_rxn,1);
for ifl = 2:FBAmodel.nt_rxn
%     ifl = 3;
    vLPmxAll{ifl} =...
    zeros(FBAmodel.nt_rxn,length(flval),FBAmodel.nt_rxn);
    vLPmnAll{ifl} = ...
    zeros(FBAmodel.nt_rxn,length(flval),FBAmodel.nt_rxn);
    Maxtarget{ifl} = zeros(1,length(flval),FBAmodel.nt_rxn);
    Mintarget{ifl} = zeros(1,length(flval),FBAmodel.nt_rxn);
    bounds.vl = zeros(FBAmodel.nt_rxn,1);    
    bounds.vl(5) = -1;
    bounds.vu = zeros(FBAmodel.nt_rxn,1);
    bounds.vu(bounds.vu==0) = 1;
    for iv = 1:length(flval)
        bounds.vl(ifl) = flval(iv);
        bounds.vu(ifl) = flval(iv);
        for imax = 7%2:FBAmodel.nt_rxn
            prxnid = imax;
            [vLPmax,vLPmin,fmax,fmin] = solveLP(FBAmodel,'P','P5',bounds,prxnid);
             if ~isempty(vLPmax)
                vLPmxAll{ifl}(:,iv,imax) = vLPmax;
                Maxtarget{ifl}(1,iv,imax) = fmax;
            end
            if ~isempty(vLPmin)
                vLPmnAll{ifl}(:,iv,imax) = vLPmin;
                Mintarget{ifl}(1,iv,imax) = fmin;
            end
            vLPmax = [];
            vLPmin = [];
            fmax = [];
            fmin = [];
        end
    end
end

%Plot Envelope
if rem(FBAmodel.nt_rxn,2) == 0
    nplots = FBAmodel.nt_rxn;
else
    nplots = FBAmodel.nt_rxn+1;
end
nrows = nplots/2;
for imax = 7%2:FBAmodel.nt_rxn
    subplot(nrows,2,imax-1);
    ylabel = sprintf('Flux %s',FBAmodel.Enzyme{imax});
    hline = plot([flval fliplr(flval)],...
                 [-Maxtarget{imax}(1,:,imax) fliplr(Mintarget{imax}(1,:,imax))],...
                 'LineWidth',2);
     set(get(gca,'YLabel'),'String',ylabel);
     set(get(gca,'XLabel'),'String','variable flux');
     axis tight;   
end

     

