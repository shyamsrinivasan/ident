function mc = readCNCfromFile(fname,model)
fileid = fopen(fname);
if fileid == -1
    fprintf('File %s cannot be opened.', fname);
    mc = [];
    return;
end

C = textscan(fileid, '%s%f',...
                     'Delimiter', '\t',...
                     'TreatAsEmpty', {'None'},...
                     'HeaderLines', 1);
mc = zeros(length(model.mets),1);
for ic = 1:length(C{1})
    tfm = strcmpi(C{1}{ic},model.mets);
    if any(tfm)
        mc(tfm) = C{2}(ic);
    end
end