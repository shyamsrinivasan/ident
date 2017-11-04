% calculate averages/std dev based on all concentrations and fluxes
% provided as a single nc*nval x 
function [avg,sigma2] = calc_avgstdev(nval,conc,nc,flux,nf,par,np)
if nargin<7
    np = 0;
end
if nargin<6
    par = [];
end
if nargin<5
    nf = 0;
end
if nargin<4
    flux = [];
end

if ~isempty(conc)
    npert = size(conc,2);
    avg_x = zeros(nc,npert);
    sigma2_x = zeros(1,npert);
    for ic = 1:nc
        avg_x(ic,:) = sum(conc(ic:nc:nc*nval,:),1)./nval;    
        sigma2_x(ic,:) = std(conc(ic:nc:nc*nval,:),0,1);
    end
else
    avg_x = [];
    sigma2_x = [];
end

if ~isempty(flux)
    npert = size(flux,2);
    avg_f = zeros(nf,npert);
    sigma2_f = zeros(1,npert);
    for ifx = 1:nf
        avg_f(ifx,:) = sum(flux(ifx:nf:nf*nval,:),1)./nval;
        sigma2_f(ifx,:) = std(flux(ifx:nf:nf*nval,:),0,1);
    end
else
    avg_f = [];
    sigma2_f = [];
end

if ~isempty(par)    
    avg_p = zeros(np,1);
    sigma2_p = zeros(np,1);
    for ip = 1:np
        avg_p(ip) = sum(par(ip,:))./nval;
        sigma2_p(ip) = std(par(ip,:));
    end
else
    avg_p = [];
    sigma2_p = [];
end

avg.avg_x = avg_x;
avg.avg_f = avg_f;
avg.avg_p = avg_p;

sigma2.sigma2_x = sigma2_x;
sigma2.sigma2_f = sigma2_f;
sigma2.sigma2_p = sigma2_p;
