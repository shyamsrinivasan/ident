% runDynamicsims
% simulate from any given/multiple saddle nodes (points on the
% bistable line)

% load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\VmaxVariation_July19.mat');

% perturb new systems from old equilibrium point to detect new steady state
% get old equilibirum point

% get original ss and continuation without perturbations
clear pvec
kEcat = 1;
KEacetate = 0.1;    % or 0.02
KFbpFBP = 0.1;
vFbpmax = 1;
Lfbp = 4e6;
KFbpPEP = 0.1;
vEXmax = 1;
KEXPEP = 0.3;
vemax = 1.1;        % for bifurcation analysis: 0.7:0.1:1.3
KeFBP = 0.45;       % or 0.45
ne = 2;             % or 2
acetate = 0.1;      % a.u acetate
d = 0.25;           % or 0.25 or 0.35
kPEPout = 0.2;
pvec = [KEacetate,KFbpFBP,Lfbp,KFbpPEP,...
        KEXPEP,vemax,KeFBP,ne,acetate,d,...
        kPEPout,kEcat,vFbpmax,vEXmax];

% systems check
givenModel = @(t,x)KotteODE(t,x,model,pvec);
fluxg = Kotte_givenFlux([M;model.PM],pvec,model);
dMdtg = givenModel(0,M);    
    
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
colorSpec = chooseColors(4,{'Green','Purple','Red','Orange'});
ac = find(strcmpi(model.mets,'ac[e]'));
npts = 1;

allxeq = zeros(length(M),npts);
allxdyn = zeros(length(M),length(tspan),npts);
allxf = zeros(length(M),npts);
allfeq = zeros(length(fluxg),npts);
allfdyn = zeros(length(fluxg),length(tspan),npts);
allpvec = pvec;
solveEquilibriumODE     

% get saddle node
[orig_saddle,orig_saddlepar] = getsaddlenode(data.s1,data.x1,5e-3);
pvec(9) = orig_saddlepar;
model.PM(ac-length(orig_saddle)) = orig_saddlepar;

% perturb saddle to get steady states - parameters fixed above
eps = 1e-4;
tspanf = 0:0.1:2000;
pival = orig_saddle+eps*[1;1;1];
[~,xeq1,~,feq1] = solveODEonly(1,pival,model,pvec,opts,tspanf);
nival = orig_saddle-eps*[1;1;1];
[~,xeq2,~,feq2] = solveODEonly(1,nival,model,pvec,opts,tspanf);
fss = [feq1 feq2];
xss = [xeq1 xeq2];

% get eigenvectors at saddle - parameters fixed above
% [~,lambda,w] = getKotteJacobian(orig_saddle,pvec,model);

% calculate separatrix from eigen vector based perturbations
% ht12fig=[];ha12=[];
% eps = 1e-4;
% tspanr = [0,-8.25];
% Line.Color = colorSpec{1};   
% Line.LineWidth = 2.0;
% % use only eig vectors with unstable directions
% iw = find(all(real(w)>0),1,'first');
% zi = orig_saddle+eps*w(:,iw);
% zj = orig_saddle-eps*w(:,iw);
% % separatrix curve
% xdynr_zi = solveODEonly(1,zi,model,pvec,opts,tspanr);
% [ht12fig,ha12] = FIGodetrajectories(real(xdynr_zi),orig_saddle,orig_saddle,2,[1 2],ht12fig,ha12,Line);
% xdynr_zj = solveODEonly(1,zj,model,pvec,opts,tspanr);
% [ht12fig,ha12] = FIGodetrajectories(real(xdynr_zj),orig_saddle,orig_saddle,2,[1 2],ht12fig,ha12,Line);

% needed variables: alliidpvec,alliidxeq,alliidfeq,tout,ap;
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
ndp = size(alliidpvec,3);
tout = tspan;


for iid = 1:ndp
    hf1 = [];
    ha1 = [];
    hf2 = [];
    ha2 = [];
    % collect points capable of mss
    if isfield(allnss,sprintf('iid%d',iid))
        msspts = find(allnss.(['iid' num2str(iid)]));
        sslps = allnss.(['iid' num2str(iid)])(msspts);
        ss = unique(sslps);
        nss = length(ss);
        allmsspts = [];
        for iss = 1:nss
            allmsspts = union(allmsspts,msspts(sslps==ss(iss)));
            ivalpts = zeros(2*nvar,npts);
            xeqpts = zeros(2*nvar,npts);
            eqid = zeros(2,npts);
            ivalid = zeros(2,npts);
            % perturbation for all points
            for ipt = 1:npts
                pvec = alliidpvec(ipt,:,iid);
                pvec(9) = orig_saddlepar;
                model.PM(ac-length(orig_saddle)) = orig_saddlepar;
                % if point not capable of mss
                if ~ismember(ipt,allmsspts)                   
                    % perturbations from ss 
                    [ivalpts,ivalid,xeqpts,eqid,hf1,ha1] = ParameterPerturbations(model,pvec,...
                        xss,ivalpts,ivalid,xeqpts,eqid,ipt,tspanf,colorSpec,opts,hf1,ha1);
                else
                    s1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).s1;
                    x1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).x1;
                    f1 =...
                    siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).f1;
                    index =...
                    cat(1,siid.(['iid' num2str(iid)]).(['pt' num2str(ipt)]).s1.index);
                    bifurcationPlot(x1,s1,f1,[4,2]);
%                     bifurcationPlot(x1,s1,f1,[4,1]);
%                     bifurcationPlot(x1,s1,f1,[4,3]); 
                    
                    % perturbations from ss 
                    [ivalpts,ivalid,xeqpts,eqid,hf1,ha1] = ParameterPerturbations(model,pvec,...
                        xss,ivalpts,ivalid,xeqpts,eqid,ipt,tspanf,colorSpec,opts,hf1,ha1);
                    
                    % get saddle node
%                     eps1 = 1e-4;
%                     saddle = [];
%                     while isempty(saddle)
%                         [saddle,saddlepar] = getsaddlenode(s1,x1,eps1);
%                         eps1 = eps1*10;
%                     end
%                     pvec(ap) = saddlepar;
%                     model.PM(ac-length(saddle)) = saddlepar;
%                     [~,lambda,w] = getKotteJacobian(saddle,pvec,model);
%                     
%                     % positive perturb around saddle
%                     ival1 = saddle+1e-2*[1;1;1];
%                     [~,xeq1] = solveODEonly(1,ival1,model,pvec,opts,tspanf);
%                     xeqpts(1:nvar,ipt) = xeq1;
%                     if xeq1(1)>xeq1(2)
%                         % if pep > fdp                        
%                         Point.MarkerFaceColor = colorSpec{1};   
%                         Point.MarkerEdgeColor = colorSpec{1}; 
%                         eqid(1,ipt) = 1;
%                     elseif xeq1(2)>xeq1(1)
%                         % if pep < fdp                        
%                         Point.MarkerFaceColor = colorSpec{2}; 
%                         Point.MarkerEdgeColor = colorSpec{2};  
%                         eqid(1,ipt) = 2;
%                     end       
%                     Point.Marker = '.';
%                     Point.MarkerSize = 20;
%                     [hf1,ha1] =...
%                     FIGmssEqIvalPerturbations(ival1,xeq1,2,[1 2],hf1,ha1,Point); 
%                     [hf2,ha2] =...
%                     FIGmssEqIvalPerturbations(ival1,xeq1,2,[1 3],hf2,ha2,Point);
                    
                    % negative pertrubation from saddle
%                     ival2 = saddle-1e-2*[1;1;1];
%                     [~,xeq2] = solveODEonly(1,ival2,model,pvec,opts,tspanf);
%                     xeqpts(nvar+1:end,ipt) = xeq2;
%                     if xeq2(1)>xeq2(2)
%                         % if pep > fdp                        
%                         Point.MarkerFaceColor = colorSpec{1};   
%                         Point.MarkerEdgeColor = colorSpec{1}; 
%                         eqid(2,ipt) = 1;
%                     elseif xeq2(2)>xeq2(1)
%                         % if pep < fdp                        
%                         Point.MarkerFaceColor = colorSpec{2}; 
%                         Point.MarkerEdgeColor = colorSpec{2};  
%                         eqid(2,ipt) = 2;
%                     end     
%                     Point.Marker = '.';
%                     Point.MarkerSize = 20;
%                     [hf1,ha1] =...
%                     FIGmssEqIvalPerturbations(ival2,xeq2,2,[1 2],hf1,ha1,Point); 
%                     [hf2,ha2] =...
%                     FIGmssEqIvalPerturbations(ival2,xeq2,2,[1 3],hf2,ha2,Point);
                end                
            end          
        end
        % plot points from xeqpts and ivalpts using ivalid and eqid after
        % normalization
        normeqpts = xeqpts./repmat(max(xeqpts,[],2),1,size(xeqpts,2));
        normival = ivalpts./repmat(max(xeqpts,[],2),1,size(ivalpts,2));
        hf1 = [];
        ha1 = [];
        plotpoints = zeros(1,size(xeqpts,2));
        Point.Marker = '.';
        Point.MarkerSize = 25;
        Point2.Marker = '.';
        Point2.MarkerSize = 25;
        for ipt = 1:size(xeqpts,2) 
            addanot.text = ['P' num2str(ipt)];
            if eqid(1,ipt)==eqid(2,ipt)
                if eqid(1,ipt)==1
                    Point.MarkerFaceColor = colorSpec{1};   
                    Point.MarkerEdgeColor = colorSpec{1}; 
                elseif eqid(1,ipt)==2
                    Point.MarkerFaceColor = colorSpec{2};   
                    Point.MarkerEdgeColor = colorSpec{2}; 
                end 
                ival1 = normival(1:nvar,ipt);
                ival2 = normeqpts(1:nvar,ipt);
                plotpoints(ipt) = 1;
            else
                if eqid(1,ipt)==1
                    Point.MarkerFaceColor = colorSpec{1};   
                    Point.MarkerEdgeColor = colorSpec{1};
                    ival1 = normival(1:nvar,ipt);
                    ival2 = normeqpts(1:nvar,ipt);
                elseif eqid(1,ipt)==2
                    Point.MarkerFaceColor = colorSpec{2};   
                    Point.MarkerEdgeColor = colorSpec{2};
                    ival1 = normival(1:nvar,ipt);
                    ival2 = normeqpts(1:nvar,ipt);
                end
                if eqid(2,ipt)==1
                    Point2.MarkerFaceColor = colorSpec{1};   
                    Point2.MarkerEdgeColor = colorSpec{1}; 
                    ival3 = normival(nvar+1:end,ipt);
                    ival4 = normeqpts(nvar+1:end,ipt);
                elseif eqid(2,ipt)==2
                    Point2.MarkerFaceColor = colorSpec{2};   
                    Point2.MarkerEdgeColor = colorSpec{2};
                    ival3 = normival(nvar+1:end,ipt);
                    ival4 = normeqpts(nvar+1:end,ipt);
                end
            end
            [hf1,ha1] =...
            FIGmssEqIvalPerturbations(ival1,ival2,2,[1 2],hf1,ha1,Point,addanot);    
            if ~plotpoints(ipt)
                [hf1,ha1] =...
                FIGmssEqIvalPerturbations(ival3,ival4,2,[1 2],hf1,ha1,Point2,addanot);    
            end          
        end
    end
end

% hobj = get(ha1,'Children');
% set(0,'CurrentFigure',hf1);
% set(hf1,'CurrentAxes',ha1);
% iter = 0;
% pid = 1;
% while 2*iter+1 <= length(hobj)
%     xdata = get(hobj(2*iter+1),'XData');
%     ydata = get(hobj(2*iter+1),'YData');
%     if rem(2*iter+1,4)==0
%         pid = pid+1;
%     end
%     ptid = [P num2str(pid)];
%     text(xdata,ydata,ptid);
%     iter = iter+1;    
% end
