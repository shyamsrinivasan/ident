sbcmp = zeros(length(model.mets),1);
prcmp = zeros(length(model.mets),1);
if any(sbid) && any(prid)
    if any(sbid([pic pie hc he h2o co2]))
            sbcmp(vmet(logical(sbid(vmet)))) = sbid(vmet(logical(sbid(vmet))));
        if any(prid([pic pie hc he h2o co2]))
            prcmp(vmet(logical(sbid(vmet)))) = prid(vmet(logical(sbid(vmet))));
            sbid([pic pie hc he h2o co2]) = 0;
            prid([pic pie hc he h2o co2]) = 0;
if model.rev(Vind(irxn))
    vflux = vfwd*mc(logical(sbcmp))^sbcmp-vbkw*mc(logical(prcmp))^prcmp