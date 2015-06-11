%Roots of the non linear integrated model 
clc
trnmodel.ext_MC = assign_extconc(batch.init{1},batch.init{2},trnmodel);
% initval = zeros(sum(ng)-ng(4),1);   
% initval(1:ng(1)) = trnmodel.SSmRNA;
% initval(ng(1)+1:ng(1)+ng(2)) = trnmodel.SSreg(1:ng(2));
% initval(ng(1)+ng(2)+1:sum(ng)-ng(4)) = trnmodel.SSreg(ng(2)+1:ng(2)+ng(3));
initval = Solution.initSS.y(:,end);
[xold] = nonlin_newtons(initval,defparval,ng,trnmodel,2000);