function printCS(cutsets,reacID)
if ~iscell(reacID)
    reacID = cellstr(reacID);
end

for i = 1:size(cutsets,1);
    fprintf('Cut set #%d',i);
    reacID(logical(cutsets(i,:)))
    fprintf('\n');
end