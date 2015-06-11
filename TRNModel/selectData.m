%function status = selectData(Solution,model,ng)
%Write select data to excel file 
function status = selectData(Solution,model,ng)
names = fieldnames(Solution);
nsim = length(names);
if ~isempty(names)
    [nvar,~] = size(Solution.(names{1}).y);    
end
%Filename
saveData.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Results\TRN Model version 4';
curr_dir = pwd;
cd(saveData.dirname);
fprintf('Writing selected data to Excel File\n');
Var = cell(nvar,1);
Var(1:ng(1)) = model.Gene;
Var(ng(1)+1:ng(1)+ng(2)) = model.Enzyme;
Var(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3)) = model.Metabolites(1:ng(3));
%Find Regulator and Regulated
[gene,reg] = find(model.RS);
reg = reg + ng(1);
indx = [gene;reg];
g_name = Var(gene);
%Add protein name
for ig = 1:length(g_name)
    var_tf = strcmpi(g_name{ig},Var);
    if any(var_tf(1:ng(1)))
        if ~any(indx==ng(1)+find(full(model.trate*var_tf(1:ng(1)))))
            indx = [indx;ng(1)+find(full(model.trate*var_tf(1:ng(1))))];
        end       
    end    
end
indx = [indx;(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3))'];
nindx = length(indx);
%Headers
wrt_cell = cell(1,nindx+1);
wrt_cell(1) = {'Time'};
wrt_cell(2:nindx+1) = Var(indx)';
for isim = 1:nsim
    sim_name = names{isim};
    npts = size(Solution.(sim_name).y(indx,:)',1);
    wrt_cell(2:npts+1,1) = num2cell(Solution.(sim_name).t/3600);
    wrt_cell(2:npts+1,2:nindx+1) = num2cell(Solution.(sim_name).y(indx,:)');   
    status = xlswrite('sim_results',wrt_cell,sim_name);
end
cd(curr_dir);
return