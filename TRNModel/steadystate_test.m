trnmodel.nvar = sum(ng)-ng(4);%%[mRNA,Protein,Metabolite]'
initval = zeros(trnmodel.nvar,1);
initval(1:ng(1)) = trnmodel.SSmRNA;%mRNA umole
% initval(1:ng(1)) = 1e-6;
initval(ng(1)+1:ng(1)+ng(2)) = trnmodel.SSreg(1:ng(2));%Protein umole
% initval(ng(1)+1:ng(1)+ng(2)) = 1e-5;
initval(ng(1)+ng(2)+1:sum(ng)-ng(4)) = trnmodel.SSreg(ng(2)+1:ng(2)+ng(3));
% initval(ng(1)+ng(2)+1:sum(ng)-ng(4)) = 1;

trnmodel.ext_MC = [1e3;10];
trnmodel.Yref = ones(trnmodel.nvar,1);

data.par = trnmodel.allpar;
data.ng = ng;%[ngene;nprot;nregp;nrecp;nmetab];
data.pdecay = 3.1e-3;
data.mdecay = 3.83e-6;
data.kcat = 1;
data.gmax = 0.2;
data.rephill = defparval.rephill;
data.ext_MC = trnmodel.ext_MC;
data.MC = initval(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3));
data.Yref = trnmodel.Yref;
data.flux = trnmodel.Vss;

initval = trnmodel.SSmRNA(1);
callsystem = @(Y)algebraicSystem(Y,data,trnmodel,FBAmodel);
options = optimset('Display','iter',...
                   'FunValCheck','on',...
                   'PlotFcns',...
                   {@optimplotfval},...
                   'TolX',1e-6);
%                    'MaxFunEvals',100000,...
%                    'MaxIter',10000,...
%                    'TolFun',1e-6,...
                   
% ,@optimplotfirstorderopt
initval = 1e-5;
% [Y,fval] = fzero(callsystem,initval,options);
% Yres = Y;
% [Y,fval] = fzero(callsystem,Yres,options);
% Yres2 = Y;
% [Y,fval] = fzero(callsystem,Yres2,options);
% Yres3 = Y;
% [Y,fval] = fzero(callsystem,Yres3,options);
% Yres4 = Y;
% [Y,fval] = fzero(callsystem,Yres4,options);
% Yres5 = Y;
% [Y,fval] = fzero(callsystem,Yres5,options);
% [Y,fval] = fsolve(callsystem,initval,options);
% Yres=Y;
% [x1,fval,exitflag,output] = fminsearch(callsystem,initval,options);

mdecay = 0.0031;
pdecay = 3.83e-6;
gr_flux = 0.8/3600;
alphac = 0.06;
alpha = 40/(233/(gr_flux*3600)^2+78);
beta = 60/(82.5/(gr_flux*3600)+148);
RNAP = 3e-5;
b = 500;
k = 1e-4;%PAR
% k = 1e-6;%NAR
c1 = ((0.06*beta/k)/(pdecay+gr_flux));
% c1 = (mdecay + gr_flux)/(RNAP*alphac);
c2 = (mdecay + gr_flux)/(RNAP);
% c2 = ((0.06*beta/k)/(pdecay+gr_flux))^2;
% c3 = alpha/(alphac*b);
% rr = c1*c2*Y^2+(c2 - alpha/b*c1 - alphac*c1)*Y - alpha/b;

Y = roots([c1*c2  (c2-alpha/b*c1-alphac*c1) -alpha/b]);
% Y = roots([c1*c2 -c2*c3 c1 -c3]);



