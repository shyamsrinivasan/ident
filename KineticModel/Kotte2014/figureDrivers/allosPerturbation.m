% build stoichioemtrc matrices
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014.txt';
cnfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\Kotte2014\Kotte2014C.txt';

% create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

% obtain conentrations from file
[mc,FBAmodel,met] = readCNCfromFile(cnfname,FBAmodel);

% run FBA
Vup_struct.ACt2r = 1;
Vup_struct.ENZ1ex = 1;
FBAmodel = FBAfluxes(FBAmodel,'fba',{'ACt2r','ENZ1ex'},Vup_struct,...
                    [find(strcmpi(FBAmodel.rxns,'FDex'))...
                     find(strcmpi(FBAmodel.rxns,'PEPex'))]);
                 
% remove metabolites held constant from consideration in the model
% integration phase
[model,pvec,newmc,cnstmet] =...
remove_eMets(FBAmodel,parameter,mc,[FBAmodel.Vind FBAmodel.Vex],...
{'enz1[c]','enz1[e]','enz[e]','ac[e]','bm[c]','bm[e]','pep[e]'});

% only initialize for varmets   
nvar = length(model.mets)-length(find(cnstmet));
M = newmc(1:nvar);
PM = newmc(nvar+1:end);
model.PM = PM;       

% call to parameter sampling script for analysis of mss
% parameters
clear pvec
k1cat = 1;
K1ac = 0.1;    % or 0.02
K3fdp = 0.1;
v3max = 1;
L3 = 4e6;
K3pep = 0.1;
v2max = 1;
K2pep = 0.3;
vemax = 1.1;        % for bifurcation analysis: 0.7:0.1:1.3
KeFDP = 0.45;       % or 0.45
ne = 2;             % or 2
acetate = 0.1;      % a.u acetate
d = 0.25;           % or 0.25 or 0.35
k4cat = 0.2;
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFDP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];

% systems check
givenModel = @(t,x)KotteODE(t,x,model,pvec);
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
tspan = 0:0.1:2000;
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
ac = find(strcmpi(model.mets,'ac[e]'));
npts = 1;
ap = 9;
allpvec = pvec;

% find equilibrium solution and run equilibrium continuation
allxeq = zeros(length(M),npts);
allxdyn = zeros(length(M),length(tspan),npts);
allxf = zeros(length(M),npts);
allfeq = zeros(length(fluxg),npts);
allfdyn = zeros(length(fluxg),length(tspan),npts);
solveEquilibriumODE

% get saddle node to get both stable steady states and get a bistable
% acetate concentration for perturbation
[orig_saddle,orig_saddlepar] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = orig_saddlepar;
model.PM(ac-length(orig_saddle)) = orig_saddlepar;

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = orig_saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = orig_saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);
xss = [xeq1 xeq2];

idp = 4;
type = 'together';
cmb = linspace(0,0.5,30)';
npts = size(cmb,1);

% set acetate conentration to saddle
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFDP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];
pvec(3) = 3e6;
% pvec(ap) = orig_saddlepar;
% model.PM(ac-length(orig_saddle)) = orig_saddlepar;

% set parameters from cmb at idp position(s)
allpvec = repmat(pvec,npts,1);
allpvec(:,idp) = cmb;

% find equilibrium point for lowest acetate 
allpvec(:,ap) = 0.01;
model.PM(ac-length(orig_saddle)) = 0.01;    
[~,allxeqlac] = solveODEonly(npts,M,model,allpvec,opts,tspan);

% and continue on acetate
[s,mssid,nss] = setupMATCONT(allxeqlac,allpvec,ap,model,fluxg,npts,1500);

%% get boundaries of acetate bistability
K3pepbsval = allpvec(mssid,idp);
acbounds = zeros(2,length(mssid)); % [min;max];
% xbounds = zeros(nvar,2*length(mssid));
mssipt = 1;
for ipt = 1:npts
    if ismember(ipt,mssid)
        index = cat(1,s.(['pt' num2str(ipt)]).s1.index);
        x1 = s.(['pt' num2str(ipt)]).x1;
%         xcont = x1(1:nvar,index);
        pcont = x1(nvar+1:end,index);
%         xbounds(:,2*mssipt-1:2*mssipt) = xcont(:,2:end-1);
        acbounds(1,mssipt) = min(pcont(:,2:end-1));
        acbounds(2,mssipt) = max(pcont(:,2:end-1));
        mssipt = mssipt+1;
    end
end
figure
plot(acbounds(1,:),K3pepbsval,'Color','b','LineWidth',2);
hold on
plot(acbounds(2,:),K3pepbsval,'Color','b','LineWidth',2);
xlabel('acetate a.u.');
ylabel('k3pep a.u.');

fprintf('Summary of k3pep Perturbation\n')
fprintf('\t\t K3pep \t acetate\t\n');
fprintf('Minimum Bistable \t %4.3f \t %4.3f\n',min(K3pepbsval),min(acbounds(1,:)));
fprintf('Maximum Bistable \t %4.3f \t %4.3f\n',max(K3pepbsval),max(acbounds(2,:)));

%% get pep and v4 vs K3pep for different acetate 
acetate = orig_saddlepar;
ap = 9;
colorSpec = chooseColors(5,{'Green','Purple','Red','Navy','HotPink'});
saddleac = zeros(npts,length(acetate));
xeqac = zeros(2*nvar,npts,length(acetate));
feqac = zeros(2*length(fluxg),npts,length(acetate));
for iac = 1:length(acetate)    
    % calculate saddle for each acetate concentration
    eps = 1e-4;
    [saddle,saddlepar,status] = eqptwrapper(s,nvar,acetate(iac),eps);
    
    % good saddle node points only
    goodsaddle = saddle(:,logical(status));
    goodsaddlepar = saddlepar(logical(status));
    saddleac(logical(status),iac) = goodsaddlepar;
    allpvec(logical(status),ap) = goodsaddlepar;
    % saddle parameter out of bifurcation bounds
    % get the one possible steady state
    oubsaddle = saddle(:,~logical(status));
    oubsaddlepar = saddlepar(~logical(status));
    saddleac(~logical(status),iac) = acetate(iac);
    allpvec(~logical(status),ap) = acetate(iac);
    
    for ipt = 1:npts
        if ismember(ipt,find(status))
            model.PM(ac-length(orig_saddle)) = saddlepar(ipt); 
            % perturb saddle to get steady states
            eps = 1e-4;                            
            pival = saddle(:,ipt)+eps*[1;1;1];
            [~,xeq1,~,feq1] =...
            solveODEonly(1,pival,model,allpvec(ipt,:),opts,tspanf);
            nival = saddle(:,ipt)-eps*[1;1;1];
            [~,xeq2,~,feq2] =...
            solveODEonly(1,nival,model,allpvec(ipt,:),opts,tspanf);            
            xeqac(nvar+1:end,ipt,iac) = xeq2;
            feqac(length(fluxg)+1:end,ipt,iac) = feq2;
        else
            % get the only possible steady state
            model.PM(ac-length(orig_saddle)) = acetate(iac);
            [~,xeq1,~,feq1] = solveODEonly(1,M,model,allpvec(ipt,:),opts,tspan);            
            xeqac(nvar+1:end,ipt,iac) = xeq1;
            feqac(length(fluxg)+1:end,ipt,iac) = feq1;
        end
        xeqac(1:nvar,ipt,iac) = xeq1;
        feqac(1:length(fluxg),ipt,iac) = feq1;
    end
    
    hc1 = figure;
    hold on
    plot(allpvec(:,idp),xeqac(1,:),'Color',colorSpec{1},'LineWidth',2);
    plot(allpvec(:,idp),xeqac(4,:),'Color',colorSpec{2},'LineWidth',2);
    xlabel('K3pep s-1');
    ylabel('PEP a.u.');
    
    hc2 = figure;
    hold on
    plot(allpvec(:,idp),feqac(5,:),'Color',colorSpec{1},'LineWidth',2);
    plot(allpvec(:,idp),feqac(10,:),'Color',colorSpec{2},'LineWidth',2);
    xlabel('K3pep s-1');
    ylabel('v4 a.u.');
end

load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\RegulationVariation_July30.mat');
% idp = [3 4] % sampled together in a mesh grid
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
ndp = size(alliidpvec,3);

if exist('tspan','var')
    tout = tspan;
end

for iid = 1:ndp
    allpvec = alliidpvec(:,:,iid);
    % unique alpha2
    alpha2_range = unique(allpvec(:,idp(1)));
    nalpha2 = length(alpha2_range);
    % max and min bistable apha1 for every alpha2
    bistable_alpha1 = zeros(2,nalpha2);
    unstable_alpha1 = zeros(2,nalpha2);
    mstable_alpha1 = zeros(2,nalpha2);
    for ialpha2 = 1:nalpha2
        startid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'first');
        endid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'last');
        if isfield(allnss,sprintf('iid%d',iid))
            % find all alpha1 that are bistable
            bistable_alpha1(1,ialpha2) =...
            max(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==4,idp(2)));
            bistable_alpha1(2,ialpha2) =...
            min(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==4,idp(2)));
            % find all alpha1 that are monostable with unstable regions
            unstable_alpha1(1,ialpha2) =...
            max(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==3,idp(2)));
            unstable_alpha1(2,ialpha2) =...
            min(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==3,idp(2)));
            % find all alpha1 that are monostable
            mstable_alpha1(1,ialpha2) =...
            max(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==0,idp(2)));
            mstable_alpha1(2,ialpha2) =...
            min(allpvec(allnss.(['iid' num2str(iid)])(startid:endid)==0,idp(2)));
        end
    end
end
% figure
% plot(log(repmat(alpha2_range',2,1)),log(mstable_alpha1),'LineStyle','none','Marker','.','MarkerSize',10);
% hold on
% plot(log(repmat(alpha2_range',2,1)),log(bistable_alpha1),'LineStyle','none','Marker','.','MarkerSize',10);
% plot(log(repmat(alpha2_range',2,1)),log(unstable_alpha1),'LineStyle','none','Marker','.','MarkerSize',10);
% xlabel('log(Lfbp)');
% ylabel('log(KFbpPEP)');

figure
plot(log(alpha2_range'),mstable_alpha1(1,:),'Color','k','LineWidth',2);
hold on
plot(log(alpha2_range'),mstable_alpha1(2,:),'Color','k','LineWidth',2);
plot(log(alpha2_range'),bistable_alpha1(1,:),'Color','r','LineWidth',2);
plot(log(alpha2_range'),bistable_alpha1(2,:),'Color','r','LineWidth',2);
plot(log(alpha2_range'),unstable_alpha1(1,:),'Color','b','LineWidth',2);
plot(log(alpha2_range'),unstable_alpha1(2,:),'Color','b','LineWidth',2);
xlabel('log(L3)');
ylabel('K3PEP');

% pick one Lfbpo point and plot bifurcation diagrams for all corresponding
% KFbpPEP values
jalpha2 = size(alpha2_range,1);
if rem(jalpha2,2)~=0
    jalpha2 = (jalpha2+1)/2;
end

% bifurcation on KfbpPEP 
ap = 4;
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
for iid = 1:ndp
    allpvec = alliidpvec(:,:,iid);
    allxeq = alliidxeq(:,:,iid);
    % unique alpha2
    alpha2_range = unique(allpvec(:,idp(1)));
    nalpha2 = length(alpha2_range);    
    for ialpha2 = 1:nalpha2
        hf1 = [];
        pvec = allpvec(ialpha2,:);
        startid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'first');
        endid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'last');
        if ialpha2 == jalpha2
            Lfbp = alpha2_range(ialpha2);
            pvec(idp(1)) = alpha2_range(ialpha2);
            alpha1_range = unique(allpvec(startid:endid,idp(2)));
            nalpha1 = length(alpha1_range);   
            ialpha1 = 0;
            while (startid+ialpha1) < endid
                pvec(idp(2)) = alpha1_range(ialpha1+1);
                xeq = allxeq(:,startid+ialpha1);
                addanot.text = ['R' num2str(ialpha1+1)];
                % continue of new parameter
                ap = idp(2);
                 % run MATCONT
                [data,y,p] = execMATCONT(xeq,pvec,ap,fluxg,model);
                if ~isempty(data) && size(data.s1,1)>2
                    hf1 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hf1,addanot);
                end
                ialpha1 = ialpha1+1;
            end       
        end              
    end
end

% bifurcation on acetate diagrams for mss values only
for iid = 1:ndp
    allpvec = alliidpvec(:,:,iid);
    % unique alpha2
    alpha2_range = unique(allpvec(:,idp(1)));
    nalpha2 = length(alpha2_range);
    % max and min bistable apha1 for every alpha2
%     bistable_alpha1 = zeros(2,nalpha2);
%     unstable_alpha1 = zeros(2,nalpha2);
%     mstable_alpha1 = zeros(2,nalpha2);
    last_endid = 0;
    for ialpha2 = 1:nalpha2
        hf1 = [];
        startid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'first');
        endid = find(allpvec(:,idp(1))==alpha2_range(ialpha2),1,'last');
        if ialpha2 == jalpha2
            Lfbp = alpha2_range(ialpha2);
            if isfield(allnss,sprintf('iid%d',iid))
                mssid = find(allnss.(['iid' num2str(iid)])(startid:endid)==4);
                mssid = mssid+last_endid;
                for imss = 1:length(mssid)
                    data = siid.(['iid' num2str(iid)]).(['pt' num2str(mssid(imss))]);
    %                 s1 = siid.(['iid' num2str(iid)]).(['pt' num2str(mssid(imss))]).s1;
    %                 x1 = siid.(['iid' num2str(iid)]).(['pt' num2str(mssid(imss))]).x1;
    %                 f1 = siid.(['iid' num2str(iid)]).(['pt' num2str(mssid(imss))]).f1;
                    hf1 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hf1);
                end
            end
        end
        last_endid = endid;        
    end
end

%% 
figure
hold on
for ipt = 9877:9925 % npts
    if ismember(ipt,mssid)
        bifpts = cat(1,s.(['pt' num2str(ipt)]).s1.index);
        for i = 1:3
            x1d = s.(['pt' num2str(ipt)]).x1(:,bifpts(i)+1:bifpts(i+1)-1);
            flux = s.(['pt' num2str(ipt)]).flux(:,bifpts(i)+1:bifpts(i+1)-1);
            eig = s.(['pt' num2str(ipt)]).f1(:,bifpts(i)+1:bifpts(i+1)-1);            
            if any(real(eig)>0)
                style = '--';
            else
                style = '-';
            end
            hl1 = line(x1d(1,:),flux(1,:),'Color',colorSpec{1},'LineStyle',style);
            hl2 = line(x1d(1,:),flux(4,:),'Color',colorSpec{2},'LineStyle',style);
            hl3 = line(x1d(1,:),flux(3,:),'Color',colorSpec{3},'LineStyle',style);
        end
        xlps = s.(['pt' num2str(ipt)]).x1(:,bifpts);
        flps = s.(['pt' num2str(ipt)]).flux(:,bifpts);
        line(xlps(1,[2,3]),...
             flps(1,[2,3]),'Color','r','Marker','.','LineStyle','none');
        line(xlps(1,[2,3]),...
             flps(4,[2,3]),'Color','r','Marker','.','LineStyle','none');
        line(xlps(1,[2,3]),...
             flps(3,[2,3]),'Color','r','Marker','.','LineStyle','none');
        legend([hl1 hl3 hl2],'v1','v3','v2');     
        drawnow
    end
end
