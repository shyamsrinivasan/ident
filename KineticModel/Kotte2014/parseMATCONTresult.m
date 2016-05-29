function [mssval,ival,fval] = parseMATCONTresult(s1,y)

if size(s1,1)>2
    xindex = cat(1,s1.index);
    xindex = xindex(2:end-1);
    iindex = xindex(1);
    findex = xindex(end);
end

mssval = y(:,xindex);
ival = y(:,iindex);
fval = y(:,findex);