% add noise to input data and calculate corresponding noisy flux
function [noisy_xss,noisy_fss] = addnoise(x,flux)
if nargin<2
    flux = [];
end
nvar = size(x,1);
if ~isempty(x)    
    nsmp = size(x,2);
end
if ~isempty(flux)
    nsmp = size(flux,2);
end
% npar = size(p,1);
% no_noise_flux = zeros(size(kotte_flux_noCAS(x,p(1,:)),1),nsmp);
% for i = 1:npar
if ~isempty(flux)
    no_noise_flux = flux; % kotte_flux_noCAS(x,p);    
else
    no_noise_flux = [];    
end
nflx = size(no_noise_flux,1);

if nsmp>1   
    pd = makedist('Normal','mu',0,'sigma',.05);
%     pd = makedist('Uniform','lower',-.05,'upper',.05);    
else
    pd = makedist('Normal','mu',0,'sigma',.05);
%     pd = makedist('Uniform','lower',-.05,'upper',.05);    
end
met_noise = random(pd,nvar,nsmp);
if ~isempty(x)
    noisy_xss = x.*(1+met_noise);
else
    noisy_xss = [];
end
% additive noise on flux
flux_noise = random(pd,nflx,nsmp);
if nargout>1 && ~isempty(no_noise_flux)
    noisy_fss = no_noise_flux.*(1+flux_noise);
else
    noisy_fss = [];
end
% noisy_fss = [];

% nonlinear noise on flux from concentration
% if npar>1
%     noisy_fss = zeros(nflx,npar);
%     for i = 1:npar
%         noisy_fss(:,i) = kotte_flux_noCAS(noisy_xss(:,i),p(i,:));
%     end
% else
%     noisy_fss = kotte_flux_noCAS(noisy_xss,p);
% end

