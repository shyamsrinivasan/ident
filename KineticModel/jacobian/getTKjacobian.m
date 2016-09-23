function allvr = getTKjacobian(model,M,irxn,flux,S,rev,pratio,thetas,thetap,...
                                drflx,fwdflx,revflx)

pie = strcmpi(model.mets,'pi[e]');

% get jacobian information
allTr = zeros(length(S(:,irxn)),1);
allDr = zeros(length(S(:,irxn)),1);
allmet = S(:,irxn)~=0;
metid = find(allmet);

% numerator dTr/dci
if rev(irxn) && all(pratio>0)
    zetarxn = fwdflx./revflx;
    for im = 1:length(metid)
        if S(metid(im),irxn)<0
            if ~isinf(zetarxn)
                allTr(metid(im)) = (zetarxn*-S(metid(im),irxn)-0)/(zetarxn-1); 
            else
                allTr(metid(im)) = 1; 
            end
        elseif S(metid(im),irxn)>0
            if ~isinf(zetarxn)
                allTr(metid(im)) = (zetarxn*0-S(metid(im),irxn))/(zetarxn-1);
            else
                allTr(metid(im)) = 1;
            end
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
    if irxn == find(strcmpi(model.rxns,'PIt2r'))
%         Kapie = 0.89; % mM
%         if metid(im) == find(pie)
%             piratio = Kapie/M(metid(im));
%             allDr(metid(im)) = allDr(metid(im))/drflx-...
%             (piratio/M(metid(im)))/(1+piratio);
%         end
    else
        allDr(metid(im)) = allDr(metid(im))/drflx;
    end
%             if vecmc(metid(im),irxn)>0
%                 allDr(metid(im)) = allDr(metid(im))/vecmc(metid(im),irxn);
%             else
%                 allDr(metid(im)) =  0;
%             end
end             

% complete differential w.r.t flux(irxn) dvr/dci
allvr = flux(irxn).*(allTr - allDr);