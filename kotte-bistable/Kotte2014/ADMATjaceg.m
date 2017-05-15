function dM = ADMATjaceg(M,Extra)

dM = Extra.fh(M);

% dM = Kotte_givenNLAE(M,model,pvec);    
% if ~isempty(model)
%     PM = cons(model.PM,M);
%     allmc = [M;PM];
% else
%     allmc = M;
% end
% d = pvec(10);
% dM = zeros(3,1);
% dM = cons(dM,M);
% flux = Kotte_givenFlux(allmc,pvec,model);
% % differential equations
% % PEP
% dM(1) = flux(1) - flux(4) - flux(5);
% % FBP
% dM(2) = flux(4) - flux(3);
% % enzymes
% % E
% dM(3) = flux(2) - d*M(3);