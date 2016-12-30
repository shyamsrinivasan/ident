function reldata = extractsaddles(s1,x1)

id = cat(1,s1.index);
if length(id)>2
    id = id(2:end-1);
end
eqcontvar = x1(end,id);
% pick smallest and largest parameter value and position
[~,minid] = min(eqcontvar);
[~,maxid] = max(eqcontvar);
minid = id(minid);
maxid = id(maxid);

% get all data between minid and maxid
reldata = x1(:,min([minid maxid]):max([minid maxid]));

