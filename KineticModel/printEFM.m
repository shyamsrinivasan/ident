function reacID = printEFM(efm,idx,ray,cnap)

reacID = cellstr(cnap.reacID(idx,:));

for i = 1:size(efm,1);
    if ray(i)
        pt_ray = 'Extreme Ray';
    else
        pt_ray = 'Extreme Point';
    end
    fprintf('EFM %d - %s',i,pt_ray);
    reacID(logical(efm(i,:)))
    fprintf('\n');
end