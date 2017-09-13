% add noise to input data and calculate corresponding noisy flux
function [noisy_xss,noisy_fss,noise] = addnoise(x,flux)
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

if ~isempty(flux)
    no_noise_flux = flux; % kotte_flux_noCAS(x,p);    
else
    no_noise_flux = [];    
end
nflx = size(no_noise_flux,1);

if nsmp>1   
    pd = makedist('Normal','mu',0,'sigma',.01);
%     pd = makedist('Uniform','lower',-.05,'upper',.05);    
else
    pd = makedist('Normal','mu',0,'sigma',.01);
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

if nargout>2
    noise = struct();
    noise.x = met_noise;
    noise.f = flux_noise;
end

