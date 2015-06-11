% function [Solution,status,indx,nfcall,totaltfinish,itime] =...
%   test4(trnmodel,defparval,ngene,tnreg,ngap,InConc,tmax,scalAbsTol,initval,tout)
%**************************************************************************
%Functionn to sole ODE model using SUNDIALS CVODES (with sensitivities)
%December 2013 - version 4.0
%December 27 2013 - version 4.1 (w/ sensitivity)
%**************************************************************************
function [Solution,status,indx,nfcall,totaltfinish,itime] =...
    integrate_sens(trnmodel,defparval,ngene,tnreg,ngap,InConc,tmax,scalAbsTol,initval,tout)
%Material balance function
% dGdt(ngene x 1) = pact(ngene x 1) - decay*G(ngene x 1);
% dPdt(ngap x 1) = translationrate(ngap x ngene)*G(ngene x 1) - decay*P(ngap x 1);
% dPmdt(tnreg-ngap x 1) = rate(tnreg-ngap x nmetab)*M(nmetab x 1) -
% decay*Pm(tnreg-ngap x 1);
if nargin < 10
    tout = 0.001;
end
if nargin < 9 || length(initval) ~= ngene+tnreg
    initval = zeros(ngene+tnreg,1);
    %initval = ones(ngene+tnreg,1);
    %initval = initval/(1E-15*6.023E+23);
end
if nargin < 8
    scalAbsTol = 1e-10;
end

data.ng(1) = ngene;
data.ng(2) = ngap;
data.ng(3) = tnreg;

data.RS = trnmodel.RS;
%data.GeneRules = trnmodel.GeneRules;
data.Protein = trnmodel.Protein;
data.Metabolite = trnmodel.Metabolite;
data.Coefficient = trnmodel.Coefficient;
data.brate = trnmodel.brate;



%data.srate is a vector in v3 of promoter activity function
%data.srate is a matrix in v4 of the promoter activity fucntion
data.srate = trnmodel.srate;

data.trate = trnmodel.trate;
data.InConc = InConc;
data.Kmax = trnmodel.Kmax;
data.Ks = trnmodel.Ks;

data.defsrate = defparval.srate;
data.decay = defparval.drate;%(in terms of Mmin-1) ver 2 = 2.32e-15 M min-1 <=> 0.14 molecules min-1  ver 1 = 0.00016 s-1
data.pdecay = defparval.pdrate;
data.defkmax = defparval.kmax;
data.defhill = defparval.rephill;

data.defparval = defparval;


Newbindaffall = {};

t0 = 0.0;
AbsTol = zeros(size(initval));
%Vector of AbsTol changed from 1e-4 to 1e-6 to 1e-8 to 1e-10
%Could provide a better result?
%Can be changed from the outside by setting to ScalAbsTol
AbsTol(AbsTol==0) = scalAbsTol;


gname = {'arcA';'fnr'};
ngnames = length(gname);
indx = zeros(length(trnmodel.Gene),1);

for igname = 1:ngnames
    indx(strcmp(gname{igname},trnmodel.Gene)) = 1;
end
% indx = [indx;zeros(tnreg)];
% T = find(indx);

totaltstart = tic;

options = CVodeSetOptions('UserData',data,...                                                    
                          'RelTol',1.e-8,...
                          'AbsTol',AbsTol,...
                          'MaxNumSteps',500);
                                   
% mondata.sol = false;
% mondata.select = find(indx);
mondata.mode = 'text';
mondata.updt = 100;
mondata.skip = 10;

options = CVodeSetOptions(options,'MonitorFn',@CVodeMonitor,...
                          'MonitorData',mondata);
                      
%No Monitor Function Yet!                      
CVodeInit(@matbalance,'Adams','Functional',t0,initval,options);

% ------------------
% FSA initialization
% ------------------

yS0 = zeros(Ny,Ns);

FSAoptions = CVodeSensSetOptions('method','Simultaneous',...
                                 'ErrControl', true,...
                                 'ParamScales', data.p(plist));

CVodeSensInit(Ns, @rhsSfn, yS0, FSAoptions);



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
    [status,t,dY] = CVode(tout,'Normal');
    stats = CVodeGetStats;
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
stats = CVodeGetStats;
CVodeFree;

totaltfinish = toc(totaltstart);
                          
%old ODE solver call
%[time, dG] = ode23(@materialbalancefun,[tstart tend],initval);
Solution.t = [Solution.t;tout];
Solution.y = [Solution.y, dY];
%Solution.Newbindaff = Newbindaffall;
% Solution.g = dY(:,1:ngene);
% Solution.p = dY(:,ngene+1:end);


function [dYdt,flag,newdata] = matbalance(t,Y,data)
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

    dYdt = zeros(data.ng(1)+data.ng(3),1);
    %tP = [P;Pm];

    for igene = 1:data.ng(1)
    
        %bindaffinity(Protein,Coefficient,RS,GeneRules,igene,protconc)
        bindaff = bindaffinity(Protein,Coefficient,RS(igene,:),...
                               GeneRules{igene},igene,[P;Pm]);
                       
        %singlepromoteractivity(Protein,srate,brate,RS,GeneRules,bindaff,igene)
%         pact = singlepromoteractivity(Protein,data.srate,data.brate,...
%                                       RS(igene,:),...
%                                       GeneRules{igene},bindaff,igene);  

        %This following function call works for v3 of the function
%         pact = singlepromoteractivity_v3(Protein,data.srate(igene),data.brate,...
%                                       RS(igene,:),...
%                                       GeneRules{igene},bindaff,igene,...
%                                       data.defhill); 
%                                   
        %In version 4 of the function data.srate is a matrix instead of a
        %vector. Hence it should be passed as data.srate(igene,:) as
        %opposed to data.srate(igene)
        pact = singlepromoteractivity_v4(Protein,data.srate(igene,:),data.brate,...
                                      RS(igene,:),...
                                      GeneRules{igene},bindaff,igene,...
                                      data.defparval);
        
                              
        dYdt(igene) = pact - decay*G(igene);
    
        %Newbindaff{igene} = bindaff;    
    end

    %Newbindaffall = [Newbindaffall,Newbindaff];
    %Solution.overallaff = [Solution.overallaff,Solution.bindaff];
    dYdt(data.ng(1)+1:data.ng(1)+data.ng(2)) = data.trate*G - decay*P;

%     prodrate = recpprod(data.InConc,nmetab,data.ng(3),...
%                         data.ng(2),data.Kmax,data.Ks);
    prodrate = recpprod(data.InConc,nmetab,data.ng(3),...
                        data.ng(2),data.defkmax,data.Ks);
                
    dYdt(data.ng(1)+data.ng(2)+1:end) = prodrate - data.pdecay*Pm;
%      dYdt(data.ng(1)+data.ng(2)+1:end) = data.defkmax*data.InConc - data.pdecay*Pm;

    flag = 0;
    newdata = [];

    function [prod] = recpprod(InConc,nmetab,tnreg,ngap,Kmax,Ks)
        %nmetab = length(trnmodel.Metabolite);   
        if nmetab == tnreg-ngap
            prod = zeros(tnreg-ngap,1);
            for iregp = 1:nmetab
%                 prod(iregp) = Kmax(ngap+iregp,iregp)*InConc(iregp);
                 prod(iregp) = Kmax*InConc(iregp);
%                 prod(iregp) = Kmax(ngap+iregp,iregp)*InConc(iregp)/...
%                               (Ks(ngap+iregp,iregp)+InConc(iregp));
            end
        end    
    end
end

end
