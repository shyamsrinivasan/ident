function allK = assignKm(alls,allp,irxn,Ksample,kfwd,krev,allK,nmodels)


if nmodels > 10
    parfor im = 1:nmodels    
        allK(im).K([find(alls);find(allp)],irxn) = Ksample(:,im);
        allK(im).kfwd(irxn) = kfwd(im);
        allK(im).krev(irxn) = krev(im);
    end
else
    for im = 1:nmodels    
        allK(im).K([find(alls);find(allp)],irxn) = Ksample(:,im);
        allK(im).kfwd(irxn) = kfwd(im);
        allK(im).krev(irxn) = krev(im);
    end
end

