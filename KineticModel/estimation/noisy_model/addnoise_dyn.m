function [noisy_xdyn,noisy_fdyn] = addnoise_dyn(xdyn,fdyn,nsmp)

noisy_xdyn = cell(nsmp,1);
noisy_fdyn = cell(nsmp,1);

for jsmp=1:nsmp
    [noisy_xdyn{jsmp},noisy_fdyn{jsmp}] = addnoise(xdyn,fdyn);
end
