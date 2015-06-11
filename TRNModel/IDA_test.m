function [Solution,status,indx,nfcall,totaltfinish,itime] =...
    IDA_test(trnmodel,defparval,ngene,tnreg,ngap,InConc,tmax,initval,tout)

%Function to be tested for DAE Systems that could be solved using IDA in
%SUNDIALS

%Worth looking at if RNA polymerase is to be included in the model!!

%Function Not tested yet!!
%December 12th 2013


%Material balance function
% dGdt(ngene x 1) = pact(ngene x 1) - decay*G(ngene x 1);
% dPdt(ngap x 1) = translationrate(ngap x ngene)*G(ngene x 1) - decay*P(ngap x 1);
% dPmdt(tnreg-ngap x 1) = rate(tnreg-ngap x nmetab)*M(nmetab x 1) -
% decay*Pm(tnreg-ngap x 1);
if nargin < 9
    tout = 0.001;
end
if nargin < 8 || length(initval) ~= ngene+tnreg
    % y0 = [1.0;0.0;0.0];
    initval = zeros(ngene+tnreg,1);
    %initval = ones(ngene+tnreg,1);
    %initval = initval/(1E-15*6.023E+23);
end

data.ng(1) = ngene;
data.ng(2) = ngap;
data.ng(3) = tnreg;

data.RS = trnmodel.RS;
data.GeneRules = trnmodel.GeneRules;
data.Protein = trnmodel.Protein;
data.Metabolite = trnmodel.Metabolite;
data.Coefficient = trnmodel.Coefficient;
data.brate = trnmodel.brate;
data.srate = trnmodel.srate;
data.trate = trnmodel.trate;
data.InConc = InConc;
data.Kmax = trnmodel.Kmax;
data.Ks = trnmodel.Ks;

data.defsrate = defparval.srate;
data.decay = defparval.drate;%(in terms of Mmin-1) ver 2 = 2.32e-15 M min-1 <=> 0.14 molecules min-1  ver 1 = 0.00016 s-1
data.pdecay = defparval.pdrate;
data.defkmax = defparval.kmax;

Newbindaffall = {};

% t0 = 0.0;
t0 = 0.0;
AbsTol = zeros(size(initval));
AbsTol(AbsTol==0) = 1.0e-4;


% gname = {'arcA';'fnr'};
% ngnames = length(gname);
% indx = zeros(length(trnmodel.Gene),1);
% 
% for igname = 1:ngnames
%     indx(strcmp(gname{igname},trnmodel.Gene)) = 1;
% end
% indx = [indx;zeros(tnreg)];
% T = find(indx);

% yp0 = [-0.04;0.04;0.0];
initdval = ??;

% 



totaltstart = tic;

% options = IDASetOptions(options,'RootsFn',@rootfn, 'NumRoots',2);
options = IDASetOptions('UserData',data,...                                                    
                          'RelTol',1.e-8,...
                          'AbsTol',AbsTol,...
                          'MaxNumSteps',500,...
                          'LinearSolver','Dense');
                      
% options = IDASetOptions('UserData', data,...
%                         'RelTol',1.e-4,...
%                         'AbsTol',[1.e-8; 1.e-14; 1.e-6],...
%                         'LinearSolver','Dense',...
%                         'JacobianFn',@djacfn);

% mondata.sol = false;
% mondata.select = find(indx);
mondata.mode = 'text';
mondata.updt = 100;
mondata.skip = 10;

% options = IDASetOptions(options,'MonitorFn',@IDAMonitor,'MonitorData',mondata);
options = IDASetOptions(options,'MonitorFn',@IDAMonitor,...
                          'MonitorData',mondata);
                      
%No Monitor Function Yet!    
% IDAInit(@resfn,t0,y0,yp0,options);
IDAInit(@matresidual,t0,initval,initdval,options);

% for i = 1:6000
% count = 0;
% countmax = 100;
flag = 1;
Solution.t = [];
Solution.y = [];
nfcall = 0;

itime = [];


while flag
    itstart = tic;
    
%     [status,t,y] = IDASolve(tout,'Normal');
    [status,t,dY] = IDASolve(tout,'Normal');
    
    stats = IDAGetStats;
    %fprintf('t = %0.2e  order = %1d step = %0.2e \n',t,stats.qlast,stats.hlast); 
    Solution.t = [Solution.t;tout];
    Solution.y = [Solution.y, dY];
    %fprintf('Solution = [%14.6e   %14.6e]',dY(T(1)),dY(T(2)));
    if status == 0 
        tout = tout*1.001;        
    end
    if tout > tmax
        flag = 0;
    end
    itfinish = toc(itstart);
    itime = [itime;itfinish];
end
stats = IDAGetStats;
IDAFree;

totaltfinish = toc(totaltstart);
                          
%old ODE solver call
%[time, dG] = ode23(@materialbalancefun,[tstart tend],initval);
Solution.t = [Solution.t;tout];
Solution.y = [Solution.y, dY];

%Model residual function for IDAS
%y - initial value y(t0)
%yp - initial value y'(t0)

% function [rr, flag, new_data] = resfn(t, y, yp, data)
% % DAE residual function
% 
% r1 = data.p(1);
% r2 = data.p(2);
% r3 = data.p(3);
% 

% 
% flag = 0;
% new_data = [];
% end

% function [rr, flag, new_data] = resfn(t, y, yp, data)
% % DAE residual function
function [rr_dYdt,flag,newdata] = matresidual(t,Y,dY,data)
    nfcall = nfcall + 1;
    decay = data.decay;
    %defsrate = data.defsrate;
    RS = data.RS;
    GeneRules = data.GeneRules;
    Protein = data.Protein;
    Coefficient = data.Coefficient;
    nmetab = length(data.Metabolite);
    %Newbindaff = cell(data.ng(1),1);

    %decay = 0.00016; %Decay Rate in s-1
    G = Y(1:data.ng(1));
    P = Y(data.ng(1)+1:data.ng(1)+data.ng(2));
    Pm = Y(data.ng(1)+data.ng(2)+1:end);
    
    dG = dY(1:data.ng(1));
    dP = dY(data.ng(1)+1:data.ng(1)+data.ng(2));
    %dPm = dY(data.ng(1)+data.ng(2)+1:end);
    

    %dYdt = zeros(data.ng(1)+data.ng(3),1);
    rr_dYdt = zeros(data.ng(1)+data.ng(3),1);
    %tP = [P;Pm];
    
    % rr(1) = -r1*y(1) + r2*y(2)*y(3) - yp(1);
    % rr(2) =  r1*y(1) - r2*y(2)*y(3) - r3*y(2)*y(2) - yp(2);
    % rr(3) = y(1) + y(2) + y(3) - 1.0;

    for igene = 1:data.ng(1)
    
        %bindaffinity(Protein,Coefficient,RS,GeneRules,igene,protconc)
        bindaff = bindaffinity(Protein,Coefficient,RS(igene,:),...
                               GeneRules{igene},igene,[P;Pm]);
                       
        %singlepromoteractivity(Protein,srate,brate,RS,GeneRules,bindaff,igene)
%         pact = singlepromoteractivity(Protein,data.srate,data.brate,...
%                                       RS(igene,:),...
%                                       GeneRules{igene},bindaff,igene);  

        pact = singlepromoteractivity_v2(Protein,data.srate(igene),data.brate,...
                                      RS(igene,:),...
                                      GeneRules{igene},bindaff,igene,...
                                      data.defsrate); 
        %Differential                      
        %dYdt(igene) = pact - decay*G(igene);
        
        %Residual
        rr_dYdt(igene) = pact - decay*G(igene) - dG(igene);
    
        %Newbindaff{igene} = bindaff;    
    end

    %Newbindaffall = [Newbindaffall,Newbindaff];
    %Solution.overallaff = [Solution.overallaff,Solution.bindaff];
    
    %Differential
    %dYdt(data.ng(1)+1:data.ng(1)+data.ng(2)) = data.trate*G - decay*P;
    
    %Residual
    rr_dYdt(data.ng(1)+1:data.ng(1)+data.ng(2)) = data.trate*G - decay*P - dP;
    
%     prodrate = recpprod(data.InConc,nmetab,data.ng(3),...
%                         data.ng(2),data.Kmax,data.Ks);
     prodrate = recpprod(data.InConc,nmetab,data.ng(3),...
                        data.ng(2),data.defkmax,data.Ks);
                
    %Differential
    %dYdt(data.ng(1)+data.ng(2)+1:end) = prodrate - data.pdecay*Pm;
    
    %Residual
    rr_dYdt(data.ng(1)+data.ng(2)+1:end) = prodrate;
%     rr_dYdt(data.ng(1)+data.ng(2)+1:end) = prodrate - data.pdecay*Pm;
    
%      dYdt(data.ng(1)+data.ng(2)+1:end) = data.defkmax*data.InConc - data.pdecay*Pm;

    flag = 0;
    newdata = [];

    function [prod] = recpprod(InConc,nmetab,tnreg,ngap,Kmax,Ks)
        %nmetab = length(trnmodel.Metabolite);   
        if nmetab == tnreg-ngap
            prod = zeros(tnreg-ngap,1);
            for iregp = 1:nmetab
%                 prod(iregp) = Kmax(ngap+iregp,iregp)*InConc(iregp);
%                  prod(iregp) = Kmax*InConc(iregp);
                prod(iregp) = Kmax(ngap+iregp,iregp)*InConc(iregp)/...
                              (Ks(ngap+iregp,iregp)+InConc(iregp));
            end
        end    
    end
end

end
