function fi = Sv(rowid,mc)                                                                                                        model,pvec,mc,rowid,rxnid)
rxnid = find(model.S(rowid,:)~=0);
flux = zeros(size(model.S,2),1);
flux(rxnid) = iflux(model,pvec,mc,flux,rxnid);

fi = model.S(rowid,:)*flux;