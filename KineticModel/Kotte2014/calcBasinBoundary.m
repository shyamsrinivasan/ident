% get intersection points with attraction basin boundary
function calcBasinBoundary(xss,lambda,w,model,pvec,opts,tspanr,tspanf,eps1)
% code similar to calc1DWus w/o forward integration from every point
% covers steps 3, 4 and 5 of algorithm
if any(real(lambda)>=0)
    w = w(:,real(lambda)>=0);
else    
    return
end

% get all xs+eps*w
neigw = size(w,2);
nss = size(xss,2);
for iw = 1:neigw
    eigw = w(:,iw);
    % xs+epsw
    for iss = 1:nss
        % check if every point on the trajectory lies within eps
        alpha1 = 1;
        eps = alpha1*eps1;
        % step 3
        bndry_pts1 = xss(:,iss) + eps*eigw;  
        % step 4 - reverse time integration from bndry_pts
        xdynr1 = solveODEonly(1,bndry_pts1,model,pvec,opts,tspanr);
        while any(all(abs(repmat(xss(:,iss),1,size(xdynr1,2))-xdynr1)>=eps))            
            if all(abs(repmat(xss(:,iss),1,size(xdynr1,2))-xdynr1)<=eps) 
                % store values of alpha and the boundary point
                pt1alpha = alpha1;
                pt1final = bndry_pts1;
            end
            alpha1 = alpha1/1.1;         
            eps = alpha1*eps;
            bndry_pts1 = xss(:,iss) + eps*eigw;             
            xdynr1 = solveODEonly(1,bndry_pts1,model,pvec,opts,tspanr);
        end
        alpha2 = 1;
        eps = alpha2*eps1;
        bndry_pts2 = xss(:,iss) - eps*eigw;
        while all((xss(:,iss)-xdynr2)>=eps)   
            % step 4 - reverse time integration from bndry_pts
            xdynr2 = solveODEonly(1,bndry_pts2,model,pvec,opts,tspanr);
            if all((xss(:,iss)-xdynr2)<=eps) 
                % store values of alpha and the boundary point
                pt2alpha = alpha2;
                pt2final = bndry_pts2;
            end
            alpha2 = alpha2/10;
            eps = alpha2*eps;
            bndry_pts2 = xss(:,iss) - eps*eigw;  
        end
        % forward intergration from the final boundary point obtained in
        % step 4
        [~,xeq1] = solveODEonly(1,pt1final,model,pvec,opts,tspanf);
    end
end

