function var = removemets(mets,rmvlst)

match = cellfun(@(x)strcmpi(rmvlst,x),mets,'UniformOutput',false);
match = cellfun(@(x)find(x),match,'UniformOutput',false);
match(cellfun(@(x)isempty(x),match)) = {0};
match = cell2mat(match);

var = mets(~logical(match));
