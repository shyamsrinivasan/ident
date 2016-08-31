% nullclines for Kotte model
function [fdp,e,allpep,allfdp,a,b,c,DX,DY,DZ] = nullclines(model,pvec)
ac = find(strcmpi(model.mets,'ac[e]'));
pvec(12) = model.PM(ac-3);

% nullclines in order of variables 1. pep 2. fdp 3. e
% fix variables whose nullclines are not calculated

% pep nullcline - vary E and fdp
options = optimoptions('fsolve','Display','iter',...
                       'TolFun',1e-10,...
                       'TolX',1e-10);
opts = optimset('Display','iter',...                       
                       'TolX',1e-10);  

e = 1;
fdp = 0:0.01:5;
x0 = 5;
allpep = zeros(length(fdp),length(e));  
for ie = 1:length(e) 
    tempe = e(ie);
    everypep = zeros(length(fdp),1);
    parfor ifd = 1:length(fdp) 
        pepfun = @(x)pep_nullcline(x,tempe,fdp(ifd),pvec);        
%         [x,fval,exitflag,output] = fsolve(pepfun,x0,options);
         [x,fval,exitflag,output] = fzero(pepfun,x0,opts);
        if exitflag > 0
            everypep(ifd) = x;        
        end
    end
    allpep(:,ie) = everypep;
end

% fdp nullcline - vary E and fdp to get pep
allfdp = zeros(length(fdp),length(e)); 
for ie = 1:length(e)    
    tempe = e(ie);
    everyfdp = zeros(length(fdp),1);
    parfor ifd = 1:length(fdp) 
        fdpfun = @(x)fdp_nullcline(x,tempe,fdp(ifd),pvec);
        [x,fval,exitflag,output] = fzero(fdpfun,x0,opts);
        if exitflag > 0
            everyfdp(ifd) = x;       
        end        
    end
    allfdp(:,ie) = everyfdp;
end

[a,b,c] = ndgrid(fdp,e,0:0.1:10);
nr = size(a,1);
nc = size(a,2);
nz = size(a,3);
a = a(:);
b = b(:);
c = c(:);
givenModel = @(x)Kotte_givenNLAE(x,model,pvec);
dx = givenModel([a';b';c']);
a = reshape(a,nr,nc,nz);
b = reshape(b,nr,nc,nz);
c = reshape(c,nr,nc,nz);
DX = reshape(dx(1,:)',nr,nc,nz);
DY = reshape(dx(2,:)',nr,nc,nz);
DZ = reshape(dx(3,:)',nr,nc,nz);
DL = sqrt((DX./5).^2+(DY./10).^2+(DZ./10).^2);
DX = DX./DL;
DY = DY./DL;
DZ = DZ./DL;            
% contour3()

% e nullcline - vary e and fdp to get pep




function dM = pep_nullcline(x,e,fdp,pvec)
dM = zeros(1,size(x,2));
flux = zeros(5,size(x,2));
flux(1,:) = pvec(1).*e.*pvec(12)./(pvec(12)+pvec(2));
flux(4,:) = pvec(7).*x./(x+pvec(8));
flux(5,:) = pvec(14)*x;
dM(1,:) = flux(1,:) - flux(4,:) - flux(5,:);

function dM = fdp_nullcline(x,e,fdp,pvec)
dM = zeros(1,size(x,2));
flux = zeros(5,size(x,2));
flux(4,:) = pvec(7).*x./(x+pvec(8));
ratio = 1+fdp./pvec(3);
flux(3,:) = pvec(4).*(ratio-1).*(ratio).^3/(ratio.^4+pvec(5).*(1+x./pvec(6)).^(-4));
dM(1,:) = flux(4,:) - flux(3,:);

function dM = e_nullcline(x,pvec)
global fdp
d = pvec(13);
dM = zeros(1,size(x,2));
flux = zeros(5,size(x,2));
flux(2,:) = pvec(9).*(1-1./(1+(pvec(10)./fdp).^pvec(11)));
dM(3,:) = flux(2,:) - d*x;