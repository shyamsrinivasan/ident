function varargout = changeparamsorder(type,pvec)

% old order
% oldpvec = [kEcat,KEacetate,KFbpFBP,vFbpmax,Lfbp,KFbpPEP,vEXmax,KEXPEP,...		
%            vemax,KeFBP,ne,acetate,d,kPEPout];
% new order
% newpvec = [KEacetate,KFbpFBP,Lfbp,KFbpPEP,KEXPEP,vemax,KeFBP,ne,acetate,...
%            d,kPEPout,kEcat,vFbpmax,vEXmax];
%
% new2oldid
new2old = [12;1;2;13;3;4;14;5;6;7;8;9;10;11];
% old2newid
old2new = [2;3;5;6;8;9;10;11;12;13;14;1;4;7];
npar = length(new2old);
switch type
    case 'new2old'
        newpvec = pvec;
        [rdim,cdim] = size(newpvec);
        if cdim ~= npar
            newpvec = newpvec';
            [rdim,cdim] = size(newpvec);
            if cdim ~= npar
                varargout{1} = [];
                return
            end
        end
        oldpvec = zeros(rdim,npar);        
        oldpvec(:,1:npar) = newpvec(:,new2old);
        varargout{1} = oldpvec;
    case 'old2new'
        oldpvec = pvec;
        [rdim,cdim,pdim] = size(oldpvec);
        newpvec = zeros(rdim,cdim,pdim);
        % check if rdim or cdim == npar        
        for ipdim = 1:pdim
            if rdim == npar 
                newpvec(1:npar,:,ipdim) = oldpvec(old2new,:,ipdim);
            elseif cdim == npar
                newpvec(:,1:npar,ipdim) = oldpvec(:,old2new,ipdim);
            else
                varargout{1} = [];
                return
            end
        end 
        newpvec = zeros(rdim,npar);
        newpvec(:,1:npar) = oldpvec(:,old2new);
%         newpvec(1) = oldpvec(2);
%         newpvec(2) = oldpvec(3);
%         newpvec(3) = oldpvec(5);
%         newpvec(4) = oldpvec(6);
%         newpvec(5) = oldpvec(8);
%         newpvec(6) = oldpvec(9);
%         newpvec(7) = oldpvec(10);
%         newpvec(8) = oldpvec(11);
%         newpvec(9) = oldpvec(12);
%         newpvec(10) = oldpvec(13);
%         newpvec(11) = oldpvec(14);
%         newpvec(12) = oldpvec(1);        
%         newpvec(13) = oldpvec(4);        
%         newpvec(14) = oldpvec(7);
        varargout{1} = newpvec;
    case 'new2oldid'
        newid = pvec;
        oldid = find(new2old(new2old==newid));
        varargout{1} = oldid;
    case 'old2newid'
        oldid = pvec;
        newid = find(old2new(old2new==oldid));
        varargout{1} = newid;
end
    
