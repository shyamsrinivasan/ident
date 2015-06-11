load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\integrated_space.mat');
% load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\Solution_intmodel.mat');
Y = Solution.initSS.y(:,end);


trnmodel.nvar = size(Y,1);
initval = zeros(trnmodel.nvar,1);
initval(end) = 0.01;%gDCW Non-zero Biomass
initval(1:ng(1)) = trnmodel.SSmRNA*initval(end);%mRNA umole
initval(ng(1)+1:ng(1)+ng(2)) = trnmodel.SSreg(1:ng(2))*initval(end);%Protein umole
initval(ng(1)+ng(2)+1:sum(ng)-ng(4)) = trnmodel.SSreg(ng(2)+1:ng(2)+ng(3));
% Y = initval;

data = struct();
data.par = trnmodel.allpar;
data.ng = ng;%[ngene;nprot;nregp;nrecp;nmetab];
data.pdecay = defparval.pdecay;
data.mdecay = defparval.mdecay;
data.rephill = defparval.rephill;
data.MC = Y(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3));
data.ext_MC = variable.MC(1:2);
data.Yref = ones(trnmodel.nvar,1);
data.kcat = 1;
data.gmax = 0.8;

%trnmodel = rmfield(trnmodel,'SSval');

obj = deriv(Y,eye(trnmodel.nvar));
y_ADMAT = integrated_ODEmodel(0,obj,data,trnmodel,FBAmodel);
Y_ADMAT = getval(y_ADMAT);
J_ADMAT = getydot(y_ADMAT);



