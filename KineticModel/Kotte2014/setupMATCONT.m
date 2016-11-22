function [s,mssid,nss,hbif1] =...
        setupMATCONT(allxeq,allpvec,ap,model,fluxg,npts,bfpts,fig,hbif1,annot)
if nargin<10
    annot = [];
end
if nargin<9
    hbif1 = [];
end
if nargin<8
    fig = 0;
end
if nargin<7
    bfpts = 800;
end

for ipt = 1:npts
    fprintf('\nIteration #%d of %d Equilibrium Continuation...\n',ipt,npts);
    
    xeq = allxeq(:,ipt);
    pvec = allpvec(ipt,:);
    
    % change to bigger number if bifurcation starts at a bigger value
%     pvec(ap) = 0.001;
    
    % run MATCONT
    [data,y,p] = execMATCONT(xeq,pvec,ap,fluxg,model,bfpts);    
    if fig && (~isempty(data) && size(data.s1,1)>2)
%         bifurcationPlot(data.flux,data.s1,data.f1,[5,3]);
%         bifurcationPlot(data.x1,data.s1,data.f1,[4,2]);
%         bifurcationPlot([data.flux;data.x1(end,:)],data.s1,data.f1,[6,5]); 
        if ~isempty(hbif1) 
            if ~isempty(annot)
                hbif1 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hbif1,annot);
            else
                hbif1 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap,hbif1);
            end
        else
            hbif1 = bifurcationPlot(data.x1,data.s1,data.f1,[4,1],ap);
        end       
    end   
    
    % save MATCONT results
    s.(['pt' num2str(ipt)]) = data;
    
    % get the mss for y and p
    if ~isempty(data)
        [yss,iyval,fyval] = parseMATCONTresult(data.s1,y);
        [pss,ipval,fpval] = parseMATCONTresult(data.s1,p);
        [fss,ifval,ffval] = parseMATCONTresult(data.s1,data.flux);
    end
        
    fprintf('Equilibrium Continuation Complete\n');
    clear pvec xeq data y p
end

% check which solutions have mss
mssid = [];
nss = zeros(npts,1);
for ipt = 1:npts
    if ~isempty(s.(['pt' num2str(ipt)]))
        s1 = s.(['pt' num2str(ipt)]).s1;
        nLP = size(s1,1)-2;
        if nLP > 0
            fprintf('Vector %d has %d Limit Points\n',ipt,nLP);
            mssid = union(mssid,ipt);
            nss(ipt) = nLP;
        end
    else
        fprintf('No convergence at %d\n',ipt);
    end
end