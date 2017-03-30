% kpepout/bifurcation diagrams
% code to generate figures showing changes in system characteristics for
% changes in v4max/kpepout
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

idp = 11;
type = 'together';
cmb = linspace(0,1,50)';
npts = size(cmb,1);

% set acetate conentration to saddle
pvec = [K1ac,K3fdp,L3,K3pep,...
        K2pep,vemax,KeFDP,ne,acetate,d,...
        k4cat,k1cat,v3max,v2max];
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
[s,mssid,nss] = setupMATCONT(@KotteMATCONT,@Kottecont_fluxcalc,allxeqlac,...
                            allpvec,ap,model,fluxg,npts,1500);

%% get boundaries of acetate bistability
k4bsval = allpvec(mssid,idp);
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
plot(acbounds(1,:),k4bsval,'Color','b','LineWidth',2);
hold on
plot(acbounds(2,:),k4bsval,'Color','b','LineWidth',2);
xlabel('acetate a.u.');
ylabel('k4cat s-1');

fprintf('Summary of k4cat Perturbation\n')
fprintf('\t\t k4cat \t acetate\t\n');
fprintf('Minimum Bistable \t %4.3f \t %4.3f\n',min(k4bsval),min(acbounds(1,:)));
fprintf('Maximum Bistable \t %4.3f \t %4.3f\n',max(k4bsval),max(acbounds(2,:)));
    

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
%         hf2 = bifurcationPlot([flux1;x1(end,:)],s1,f1,[6,5],ap,hf2,addanot);
%         hf2 = bifurcationPlot(x1,s1,f1,[1,4],ap,hf2,addanot);
    end
end

%% bifurcation plot for fluxes w.r.t acetate
for ipt = 1:npts
    
end

%% get pep and v4 vs k4cat for different acetate 
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
    plot(allpvec(:,idp),xeqac(1,:),'Color',colorSpec{1},'LineWidth',2,'Marker','s','MarkerSize',11);
    plot(allpvec(:,idp),xeqac(4,:),'Color',colorSpec{2},'LineWidth',2,'Marker','o','MarkerSize',8);
    xlabel('k4cat s-1');
    ylabel('PEP a.u.');
    
    hc2 = figure;
    hold on
    plot(allpvec(:,idp),feqac(5,:),'Color',colorSpec{1},'LineWidth',2,'Marker','s','MarkerSize',11);
    plot(allpvec(:,idp),feqac(10,:),'Color',colorSpec{2},'LineWidth',2,'Marker','o','MarkerSize',8);
    xlabel('k4cat s-1');
    ylabel('v4 a.u.');
end
%% continue on acetate for given kcat to get all stable states
pvec = allpvec(18,:);
ap = 9;
pvec(ap) = 0.01;
model.PM(ac-length(orig_saddle)) = 0.01;
% pvec(11) = 0.2;
iguess = [0.1 0.1 0.1];
[~,xeq1,~,feq1] = solveODEonly(1,iguess',model,pvec,opts,tspanf);
lambda = linspace(0.01,2.5,500);
% options = optimoptions(@fsolve,'TolX',1e-6,'TolFun',1e-6,'Display','iter');
[contpt,contflux,conteig] = continuation(model,pvec,lambda,ap,xeq1,opts);
%% plot data form previous cell
figure
hold on
line(contpt(1,:),contflux(1,:),'LineStyle','none','Marker','.','Color',colorSpec{1});
line(contpt(1,:),contflux(4,:),'LineStyle','none','Marker','.','Color',colorSpec{2});
line(contpt(1,:),contflux(5,:),'LineStyle','none','Marker','.','Color',colorSpec{3});

%% get saddle close to original acetate for all bistable systems parameters
[newsaddle,saddleparpts,saddlestatus,saddleid] =...
geteqcontpt(s,orig_saddlepar,1e-3);
saddleflux = zeros(5,npts);
for ipt = 1:npts
    if ismember(ipt,mssid)
        pvec(ap) = saddleparpts(ipt);
        model.PM(ac-length(orig_saddle)) = saddleparpts(ipt);
        pvec(11) = allpvec(ipt,11);
        saddleflux(:,ipt) =...
        Kotte_givenFlux([newsaddle(:,ipt);model.PM],pvec,model);
    end
end

%% plot all saddles v4 vs pep
figure
hold on
for ipt = 1:npts
    if ismember(ipt,mssid)
       line(newsaddle(1,ipt),saddleflux(4,ipt),'Marker','.','Color','r');
       drawnow
    end
end

%% % plot all available steady state pep concentrations for a 
% given parameter V4max against all available fluxes
% choose ipt = 1, 5, 11, 15 and 17 (all points 1:17 are bistable)
% figure
hold on
for ipt = 17:17 % npts
    if ismember(ipt,mssid)
        bifpts = cat(1,s.(['pt' num2str(ipt)]).s1.index);
        for i = 1:3
            x1 = s.(['pt' num2str(ipt)]).x1(:,bifpts(i)+1:bifpts(i+1)-1);
            flux = s.(['pt' num2str(ipt)]).flux(:,bifpts(i)+1:bifpts(i+1)-1);
            eig = s.(['pt' num2str(ipt)]).f1(:,bifpts(i)+1:bifpts(i+1)-1);            
            if any(real(eig)>0)
                style = '--';
            else
                style = '-';
            end
            hl1 = line(x1(1,:),flux(1,:),'Color',colorSpec{1},'LineStyle',style);
            hl2 = line(x1(1,:),flux(4,:),'Color',colorSpec{2},'LineStyle',style);
            hl3 = line(x1(1,:),flux(5,:),'Color',colorSpec{3},'LineStyle',style);
        end
        if saddlestatus(ipt)
            line(newsaddle(1,ipt),saddleflux(1,ipt),'Color','r','Marker','.');
            line(newsaddle(1,ipt),saddleflux(4,ipt),'Color','r','Marker','.');
            line(newsaddle(1,ipt),saddleflux(5,ipt),'Color','r','Marker','.');
            drawnow
        end
        xlps = s.(['pt' num2str(ipt)]).x1(:,bifpts);
        flps = s.(['pt' num2str(ipt)]).flux(:,bifpts);
        line(xlps(1,[2,3]),...
             flps(1,[2,3]),'Color','k','Marker','.','LineStyle','none');
        line(xlps(1,[2,3]),...
             flps(4,[2,3]),'Color','k','Marker','.','LineStyle','none');
        line(xlps(1,[2,3]),...
             flps(5,[2,3]),'Color','k','Marker','.','LineStyle','none');
        legend([hl1 hl3 hl2],'v1','v4','v2');     
        drawnow
    end
end

%% % plot all available steady state pep concentrations for a 
% given parameter V4max against all available concentrations
% choose ipt = 1, 5, 11, 15 and 17
% figure
hold on
for ipt = 5:5 % npts
    if ismember(ipt,mssid)
        bifpts = cat(1,s.(['pt' num2str(ipt)]).s1.index);
        for i = 1:3
            x1 = s.(['pt' num2str(ipt)]).x1(:,bifpts(i)+1:bifpts(i+1)-1);
            flux = s.(['pt' num2str(ipt)]).flux(:,bifpts(i)+1:bifpts(i+1)-1);
            eig = s.(['pt' num2str(ipt)]).f1(:,bifpts(i)+1:bifpts(i+1)-1);            
            if any(real(eig)>0)
                style = '--';
            else
                style = '-';
            end
            hl1 = line(x1(1,:),x1(2,:),'Color',colorSpec{1},'LineStyle',style);
            hl2 = line(x1(1,:),x1(3,:),'Color',colorSpec{2},'LineStyle',style);           
        end
        xlps = s.(['pt' num2str(ipt)]).x1(:,bifpts);
        flps = s.(['pt' num2str(ipt)]).flux(:,bifpts);
        line(xlps(1,[2,3]),...
             xlps(2,[2,3]),'Color','k','Marker','.','LineStyle','none');
        line(xlps(1,[2,3]),...
             xlps(3,[2,3]),'Color','k','Marker','.','LineStyle','none');        
        legend([hl1 hl2],'fdp','E');     
        drawnow
    end
end
