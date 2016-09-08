function allvr = getCKjacobian(irxn,S,SI,sratio,pratio,allac,allib,...
                               thetas,thetap,fwdflx,revflx,nrflx,drflx)
% integrate jacobian calculation with CKinetics.m
% call from within CKinetics to calculate jacobians using information from
% CKinetics
% not complete September 2016

% jacobian information
allTr = zeros(length(S(:,irxn)),1);
allDr = zeros(length(S(:,irxn)),1);
allfr = zeros(length(SI(:,irxn)),1);
allDreg = zeros(length(SI(:,irxn)),1);
% allvr = zeros(length(S(:,irxn)),1);
allmet = S(:,irxn)~=0;
metid = find(allmet);

% numerator dTr/dci
if rev(irxn) && all(pratio>0)
    zetarxn = fwdflx./revflx;
    for im = 1:length(metid)
        if S(metid(im),irxn)<0
            allTr(metid(im)) = (zetarxn*-S(metid(im),irxn)-0)/(zetarxn-1);                    
        elseif S(metid(im),irxn)>0
            allTr(metid(im)) = (zetarxn*0-S(metid(im),irxn))/(zetarxn-1);                    
        end
%                 if vecmc(metid(im),irxn)>0
%                     allTr(metid(im)) = allTr(metid(im))/vecmc(metid(im),irxn);
%                 else
%                     allTr(metid(im)) = 0;
%                 end
    end
elseif ~rev(irxn)
    for im = 1:length(metid)
        if S(metid(im),irxn)<0
            allTr(metid(im)) = -S(metid(im),irxn);
        end
%                 if vecmc(metid(im),irxn)>0
%                     allTr(metid(im)) = allTr(metid(im))/vecmc(metid(im),irxn);
%                 else
%                     allTr(metid(im)) = 0;
%                 end
    end
end

% denominator dDr/dci        
for im = 1:length(metid)
    if S(metid(im),irxn)<0                
        allDr(metid(im)) = thetas*-S(metid(im),irxn)+0*thetap;
    elseif S(metid(im),irxn)>0                
        allDr(metid(im)) = thetas*0+thetap*S(metid(im),irxn);
    end      
%             if vecmc(metid(im),irxn)>0
%                 allDr(metid(im)) = allDr(metid(im))/vecmc(metid(im),irxn);
%             else
%                 allDr(metid(im)) = 0;
%             end
end

% denominator specific regulation dDreg/dci
% specific activation   
if any(SI(:,irxn)>0)
%             allDreg(logical(spac)) = spac(logical(spac)).*...
%                                      (1./spiratio).^spib(logical(spac))./...
%                                      vecmc(logical(spac),irxn);
end

% specific inhibition
if any(SI(:,irxn)<0)
%             allDreg(logical(spib)) = spib(logical(spib)).*...
%                                      (spiratio).^spib(logical(spib))./...
%                                      vecmc(logical(spib),irxn);
end

allDr = allDr./drflx;
allDreg = allDreg./drflx;

% partial/essential activation dacflx/dci        
if any(SI(:,irxn)>0)
    alphaA = 1./(1+acratio);
    betaA = acratio.*alphaA;
    allfr(logical(allac)) =...
    (1-rhoA)./(1-rhoA+rhoA.*betaA.^(-allac(logical(allac)))).*...
    allac(logical(allac)).*alphaA;            
end

% partial/essential inhibition dibflx/dci        
if any(SI(:,irxn)<0)
    alphaI = 1./(1+ibratio);
    betaI = ibratio.*alphaI;
    allfr(logical(allib)) = allfr(logical(allib))...
    -(1-rhoI)./(1-rhoI+rhoI.*alphaI.^(-allib(logical(allib)))).*...
    allib(logical(allib)).*betaI;
end

% collect common ids between allac and allid
%         if any(SI(:,irxn)>0) && any(SI(:,irxn)<0)
%             cmnreg = union(find(allac),find(allib));
%             % remove common ids from allac and allib
%             allac = allac(~ismember(allac,cmnreg));
%             allib = allib(~ismember(allib,cmnreg));
%             allfr(cmnreg) = allfr(cmnreg)./vecmc(cmnreg,irxn);
% %             allfr(vecmc(cmnreg,irxn)<=0) = 0;
%         end
%         if any(SI(:,irxn)>0)
%             if ~isempty(allac)
%                 allfr(logical(allac)) = allfr(logical(allac))./vecmc(logical(allac),irxn);
% %                 allfr(
%             end
%         end
%         if any(SI(:,irxn)<0)
%             if ~isempty(allib)
%                 allfr(logical(allib)) = allfr(logical(allib))./vecmc(logical(allib),irxn);
%             end  
%         end

% complete differential w.r.t flux(irxn) dvr/dci
allvr = flux(irxn).*(allfr + allTr - allDr - allDreg);
DVX(:,irxn) = allvr;        