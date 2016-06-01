% function [mc,model] = readCNCfromFile(fname,model)
% read concentrations from text file
function [mc,model,met] = readCNCfromFile(fname,model)
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
    tfm = strcmpi(strtrim(C{1}{ic}),model.mets);
    if any(tfm)
%         display(C{1}(ic));
        mc(tfm) = C{2}(ic);
        if isfield(model,'MClow')
            model.MClow(tfm) = C{2}(ic);
        end
        if isfield(model,'MChigh')
            model.MChigh(tfm) = C{2}(ic);
        end
    end
end

%get metabolite structure suitable for input to iconcetration.m
vmet = {'h2o[c]','h2o[e]','o2[e]','o2[c]','pi[e]','h[e]','h[c]','pyr[e]'};%,
mets_idx = cellfun(@(x)strcmpi(x,vmet),model.mets,'UniformOutput',false);
vmet_id = cellfun(@(x) find(x),mets_idx,'UniformOutput',false);
% vmet_id = cell2mat(vmet_id);
vmet_id(cell2mat(cellfun(@(x)isempty(x),vmet_id,'UniformOutput',false))) = {0};
vmet_id = cell2mat(vmet_id);
vmet = cellfun(@(x)strrep(x,'[c]','_c'),vmet,'UniformOutput',false);
vmet = cellfun(@(x)strrep(x,'[e]','_e'),vmet,'UniformOutput',false);
met = struct();
for im = 1:length(vmet_id)
    if any(vmet_id(im))
%         display(vmet{vmet_id(im)});
%         display('\n');
        met.(vmet{vmet_id(im)}) = mc(ismember(vmet_id,vmet_id(im)));
    end
end
    

