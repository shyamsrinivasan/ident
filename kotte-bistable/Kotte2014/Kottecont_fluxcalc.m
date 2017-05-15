function flux = Kottecont_fluxcalc(xeq,x1,p,fluxg,model,pvec,ap)

ac = find(strcmpi(model.mets,'ac[e]'));
flux = zeros(length(fluxg),size(x1,2));
if ~isempty(x1)
    for icp = 1:size(x1,2)
        pvec(ap) = p(icp);
        model.PM(ac-length(xeq)) = p(icp);
        flux(:,icp) = Kotte_givenFlux([x1(1:length(xeq),icp);model.PM],pvec,model);
    end
end
