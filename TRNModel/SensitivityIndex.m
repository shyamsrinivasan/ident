function [S,ST] = SensitivityIndex(yA,yB,yC,Ns)
S = zeros(sum(ng)+1,Ns);
ST = zeros(sum(ng)+1,Ns);
for j = 1:Ns
    f02 = sum(yA,2);
    nr = diag(yA*yC(:,:,j)') - f02;
    dr = diag(yA*yA') - f02;
    nrt = diag(yB*yC(:,:,j)') - f02;    
    S(:,j) = nr./dr;
    ST(:,j) = nrt./dr;
end
return