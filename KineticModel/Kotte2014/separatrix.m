% nullclines for Kotte model
% load('C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KotteModel\VmaxVariation_Jun01.mat');
npts = size(alliidpvec,1);
nvar = size(alliidxeq,1);
ndp = 1;
colorSpec = chooseColors(4,{'Green','Purple','Red','Orange'});
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
ac = find(strcmpi(model.mets,'ac[e]'));


for idp = 1:ndp   
    
    if isfield(allnss,sprintf('iid%d',idp));
        msspts = find(allnss.(['iid' num2str(idp)]));
        sslps = allnss.(['iid' num2str(idp)])(msspts);
        ss = unique(sslps);
        nss = length(ss);
        
        for iss = 1:nss
            allmsspts = msspts(sslps==ss(iss));
            allisspvec = alliidpvec(allmsspts,:,idp);
            allissxeq = alliidxeq(:,allmsspts,idp);
            nmsspts = length(allmsspts);
            
            % collect corresponding limit points
            for ipt = 1:nmsspts
                index =...
                cat(1,siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).s1.index);
                s1 =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).s1;
                x1ind =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).x1(:,index);
                x1 =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).x1;
                flux1 =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).flux;
                eigind =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).f1(:,index);
                f1 =...
                siid.(['iid' num2str(idp)]).(['pt' num2str(allmsspts(ipt))]).f1;
                xLPval = x1ind(1:nvar,:);
                pLPval = x1ind(nvar+1:end,:);
                pvec = allisspvec(1,:);
                
                % find the 2 or more steady states from the LPs
                LPxeq = [];
                for it = 1:size(xLPval,2)
                    ival = xLPval(:,it); 
                    pvec(ap) = pLPval(it);
%                     model.PM(ac-length(ival)) = pLPval(it);
                    [~,xeq] =...
                    solveODEonly(1,ival,model,pvec,opts,tout); 
                    if ~isempty(LPxeq)
                        if ~any(abs(LPxeq-repmat(xeq,1,size(LPxeq,2)))<=1e-8)
                            LPxeq = [LPxeq xeq];
                        end
                    else
                        LPxeq = xeq;
                    end                
                end
                acetate = 0.06:0.01:0.6;
                [x,y,z] = ndgrid(0:0.2:2.5,0:0.3:5,0:0.1:3.5);
                nr = size(x,1);
                nc = size(x,2);
                nz = size(x,3);
                x = x(:);
                y = y(:);
                z = z(:);
%                 for i = 1:length(acetate)
                    pvec(ap) = 0.1;
                    model.PM(ac-length(ival)) = 0.1;
                    givenModel = @(x)Kotte_givenNLAE(x,model,pvec);
                    dx = givenModel([x';y';z']);
                    x = reshape(x,nr,nc,nz);
                    y = reshape(y,nr,nc,nz);
                    z = reshape(z,nr,nc,nz);
%                     DX = dx(1,:)';
%                     DY = dx(2,:)';
%                     DZ = dx(3,:)';
                    DX = reshape(dx(1,:)',nr,nc,nz);
                    DY = reshape(dx(2,:)',nr,nc,nz);
                    DZ = reshape(dx(3,:)',nr,nc,nz);
                    DL = sqrt((DX./2.5).^2+(DY./5).^2+(DZ./3.5).^2);
                    DX = DX./DL;
                    DY = DY./DL;
                    DZ = DZ./DL;
                    for iz = 1:size(DX,3)
                        contour3(reshape(y(:,:,iz),186,131),reshape(x(:,:,iz),186,131),reshape(DX(:,:,iz),186,131));
                        hold on
                    end
%                 end
            end
        end
    end
end
