function noisy_ss = addnoise_wrapper(pt_sol_id,solstruct)

nptsol = length(pt_sol_id);
noisy_ss = struct([]);
for j = 1:nptsol
    xss = solstruct(pt_sol_id(j)).xss;
    pt_p = solstruct(pt_sol_id(j)).odep;    
    [noisy_xss,noisy_fss] = addnoise(xss,pt_p);
    noisy_ss(j).xss = noisy_xss;
    noisy_ss(j).fss = noisy_fss;
end