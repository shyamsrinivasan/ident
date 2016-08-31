% parse results from MATCONT using onformation from the s1 sructure
% mss - one or more LP information
% ival - initial value of y at y(0), point at which continuation is tarted
% fval - final value of y at y(Inf), point at which continuation ends
function [mssval,ival,fval] = parseMATCONTresult(s1,y)

xindex = cat(1,s1.index);
if size(s1,1)>2    
    mssindex = xindex(2:end-1);
    mssval = y(:,mssindex);
else
    mssval = [];
end
iindex = xindex(1);
findex = xindex(end);

ival = y(:,iindex);
fval = y(:,findex);