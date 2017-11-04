% vemax bifurcation diagrams and perturbation results for design
runKotte

% get saddle node to get both stable steady states and get a bistable
% acetate concentration for perturbation
[orig_saddle,orig_saddlepar,saddleid] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(ap) = orig_saddlepar;
model.PM(ac-length(orig_saddle)) = orig_saddlepar;
saddleflux = Kotte_givenFlux([orig_saddle;model.PM],pvec,model);

% perturb saddle to get steady states
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = orig_saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = orig_saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);
xss = [xeq1 xeq2];

idp = 6;
type = 'together';
cmb = linspace(0.1,30,30)';
npts = size(cmb,1);

% set acetate conentration to saddle
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFDP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];
    
% set parameters from cmb at idp position(s)
allpvec = repmat(pvec,npts,1);
allpvec(:,idp) = cmb;

% find equilibrium point for lowest acetate 
allpvec(:,ap) = 0.001;
model.PM(ac-length(orig_saddle)) = 0.001;    
[~,allxeqlac] = solveODEonly(npts,M,model,allpvec,opts,tspan);

% and continue on acetate
[s,mssid,nss] = setupMATCONT(@KotteMATCONT,@Kottecont_fluxcalc,allxeqlac,allpvec,ap,model,fluxg,npts,6000);

%% get boundaries of acetate bistability
vemax = allpvec(mssid,idp);
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
plot(acbounds(1,:),vemax,'Color','b','LineWidth',2);
hold on
plot(acbounds(2,:),vemax,'Color','b','LineWidth',2);
xlabel('acetate a.u.');
ylabel('vemax s-1');

fprintf('Summary of vemax Perturbation\n')
fprintf('\t\t vemax \t acetate\t\n');
fprintf('Minimum Bistable \t %4.3f \t %4.3f\n',min(vemax),min(acbounds(1,:)));
fprintf('Maximum Bistable \t %4.3f \t %4.3f\n',max(vemax),max(acbounds(2,:)));

%% plot bifurcation (w/ 2 LPs) on acetate
hf1 = [];
hf2 = [];
for ipt = 1:npts
    addanot.text = ['R' num2str(ipt)];
    if ismember(ipt,mssid)
        s1 = s.(['pt' num2str(ipt)]).s1;
        x1 = s.(['pt' num2str(ipt)]).x1;
        f1 = s.(['pt' num2str(ipt)]).f1;
        flux1 = s.(['pt' num2str(ipt)]).flux;
        index = cat(1,s.(['pt' num2str(ipt)]).s1.index);
        hf1 = bifurcationPlot(x1,s1,f1,[4,1],ap,hf1,addanot);  
    end
end

%% get pep and v4 vs vemax for different acetate 
acetate = orig_saddlepar; % orig_saddlepar;
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
    xlabel('vemax s-1');
    ylabel('PEP a.u.');
    
    hc2 = figure;
    hold on
    plot(allpvec(:,idp),feqac(5,:),'Color',colorSpec{1},'LineWidth',2);
    plot(allpvec(:,idp),feqac(10,:),'Color',colorSpec{2},'LineWidth',2);
    xlabel('vemax s-1');
    ylabel('v4 a.u.');
end



    