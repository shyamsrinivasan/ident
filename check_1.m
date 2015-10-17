function check_1(model,x)
A = model.A(:,1:length(model.mets));
x = x(1:length(model.mets));
Gamma = A*x;

for i=1:length(model.rxns)    
%     if abs(Gamma(i)-log(model.Keq(i)))<1e-7
        fprintf('%d. Reaction %s Gamma/Keq %3.4g Vss %3.4g\n',i,model.rxns{i},...
                                               Gamma(i)-log(model.Keq(i)),...
                                               model.Vss(i));
%     end        
end