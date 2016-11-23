% function 
% wrapper around geteqcontpt for array of points
function [allsaddle,allsaddlepar,allstatus] = eqptwrapper(s,nvar,acetate,eps)
npts = length(fieldnames(s));
allsaddle = zeros(nvar,npts);
allsaddlepar = zeros(1,npts);
allstatus = zeros(1,npts);
for ipt = 1:npts
%     epsset = eps;
%     saddle = [];
%     saddlepar = [];
%     status = -1;
%     while status<0
        s1.pt1 = s.(['pt' num2str(ipt)]);
        [saddle,saddlepar,status] = geteqcontpt(s1,acetate,eps);
%         epsset = epsset*10;
%     end
    if status        
        allsaddle(:,ipt) = saddle;
        allsaddlepar(ipt) = saddlepar;        
    end
    allstatus(ipt) = status;
end
