%function status = selectData(Solution,model,ng)
%Write select flux data to excel file 
function status = selectFlux_kmodel(Solution,model,sheet)
[nvar,npts] = size(Solution.y);    
% nint_metab = model.nint_metab;
nt_rxn = model.nt_rxn;
%Filename
saveData.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Results\Kmodel';
curr_dir = pwd;
cd(saveData.dirname);
fprintf('Writing selected data to Excel File\n');
%Column Headers
wrt_cell = cell(1,nvar+2);
wrt_cell(2) = {'Initial SS'};
ipertb = 1;
while ipertb <= nvar
    wrt_cell(ipertb+2) = {sprintf('Perturbation #%d',ipertb)};    
    ipertb = ipertb+1;
end
%Row Headers
wrt_cell(2:npts+1,1) = model.Enzyme(1:nt_rxn);
%Data
wrt_cell(2:npts+1,2) = num2cell(Solution.t);
wrt_cell(2:npts+1,3:nvar+2) = num2cell(Solution.y');
% wrt_cell(2:npts+1,1) = num2cell(Solution.t/3600);
% wrt_cell(2:npts+1,2:nvar+1) = 
status = xlswrite('sim_results',wrt_cell,sheet);
cd(curr_dir);
return