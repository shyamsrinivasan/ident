function [pts,assignFlag,delGr,vCorrectFlag] = ACHRmetSampling(model,nFiles,nptsPerFile,stepsPerPnt)
if nargin<4
    stepsPerPnt = 199;
end

if nargin<3
    nptsPerFile=1000;
end

if nargin<2
    nFiles = 10;
end

%distance closest constraint
maxMinTol = 1e-9;

%ignore directions tolerance
uTol = 1e-9;

%directions too close to boundary
dTol = 1e-14;

%create warmup points for ACHR - temporary call
bounds = setupMetLP(model);
if size(bounds.A,2)==length(bounds.mets)
    %setup slack problem
    bounds = setupSlackVariables(bounds);
    warmUpPts = createWarmupPoints(model,bounds,2000);   
end

lb = separate_slack(bounds.lb,model,bounds);
ub = separate_slack(bounds.ub,model,bounds);

[nmets,npts] = size(warmUpPts);

%centre point for the space
centre = mean(warmUpPts,2);

%set centre point as start point
prevPt = centre;

totalStepcnt = 0;

%nFiles - # files created to store samples
for i=1:nFiles
    
    pts = zeros(nmets,nptsPerFile);
    
    ptcnt = 1;
    while ptcnt<=nptsPerFile
        
        %random step size
        randStepSize = rand(stepsPerPnt,1);
        
        stepcnt = 1;
        while stepcnt<=stepsPerPnt
            
            %pick random warmup point
            randPntID = ceil(npts*rand);
            randPnt = warmUpPts(:,randPntID);
            
            %direction from centre point to warmup point
            u = randPnt-centre;
            u = u/norm(u);
            
            %distance to upper and lower bound
            distUb = ub-prevPt;
            distLb = prevPt-lb;
            
            %determine if too close to boundary
            validDir = ((distUb>dTol) & (distLb>dTol));
            
            if ~any(validDir)
                continue
            end
            
            %find positive and negative directions
            posDir = find(u(validDir)>uTol);
            negDir = find(u(validDir)<-uTol);
            
            %find all max and min step sizes
            maxStepTemp = distUb(validDir)./u(validDir);
            minStepTemp = -distLb(validDir)./u(validDir);
            
            maxStepVec = [maxStepTemp(posDir);minStepTemp(negDir)];
            minStepVec = [minStepTemp(posDir);maxStepTemp(negDir)];
            
            %true max and min step sizes
            maxStep = min(maxStepVec);
            minStep = max(minStepVec);
            
            %Find new direction if getting too close to boundary
            if (abs(minStep) < maxMinTol & abs(maxStep) < maxMinTol) | (minStep > maxStep)
                fprintf('Warning\n');
                continue
            end
            
            %obtain random step distance
            stepDist = randStepSize(stepcnt)*(maxStep-minStep)+minStep;
            
            %advance to next point
            curPt = prevPt + stepDist*u;
            
            %Reproject current point and go to the next step
%             if mode(totalStepcnt,10)==0
%                 if full(max(max(abs(model.S*curPt))))>1e-9
%                     curPt = N*(N'*curPt);
%                 end
%             end
            
            %print errors to file
            
            %print step information
%             overInd = find(ub-curPt>0);
%             underInd = find(curPt-lb<0);
            
            if (any((ub-curPt) < 0) || any((curPt-lb)< 0))
               curPt(ub-curPt>0) = ub(ub-curPt>0);
               curPt(curPt-lb<0) = lb(curPt-lb<0);             
            end
            
%             if mod(totalStepcnt,2000)==0
%                 %write to file
%             end
            
            prevPt = curPt;
            stepcnt = stepcnt+1;
            
            %count toal number of steps
            totalStepcnt = totalStepcnt+1;
            
            %recalculate centre point
            centre = ((npts+totalStepcnt)*centre+curPt)/(npts+totalStepcnt+1);
            
        end %steps per point
        
        %add current point to points
        pts(:,ptcnt) = curPt;
        ptcnt = ptcnt+1;
    end %points per cycle
    
    %re-assign concentrations to model.mets
    [pts,assignFlag,delGr,vCorrectFlag] = assignConc(pts,model,bounds); 
    [delGr,assignFlux] = assignRxns(delGr,model,bounds);
    
    mc = exp(pts);
    mc(pts==0)=0;
    pts = mc;
    
    %re-checking delGr and concentrations for thermodynamic feasibility
    vSpl = zeros(length(model.rxns),1);
    for ir = 1:length(bounds.rxns)
        vSpl(ir) = length(find(vCorrectFlag(ir,:)));
        if all(vCorrectFlag(ir,:))
            fprintf('%d. %s\n',ir,model.rxns{ir});
        end
    end
    for is = 1:length(pts(1,:))
        if all(vCorrectFlag(:,is))
            fprintf('Sample %d\n',is);
        end
    end           
          
%     [delGr
%     for ipt = 1:length(pts(1,:))
%         if ~isempty(pts(:,ipt))
%             
%         end
%     end
    %save current points to file
end
            