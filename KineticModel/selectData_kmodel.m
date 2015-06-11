%function status = selectData(Solution,model,ng)
%Write select concentration data to excel file 
function status = selectData_kmodel(Solution,model,sheet)
[nvar,npts] = size(Solution.y);    
nint_metab = model.nint_metab;
%Filename
saveData.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Results\Kmodel';
curr_dir = pwd;
cd(saveData.dirname);
fprintf('Writing selected data to Excel File\n');
%Headers
wrt_cell = cell(1,nvar+1);
wrt_cell(1) = {'Time'};
ipertb = 0;
while ipertb < (nvar/nint_metab)
    wrt_cell(nint_metab*ipertb+2:nint_metab*(ipertb+1)+1) =...
    model.Metabolites(1:nint_metab);
    ipertb = ipertb+1;
end
wrt_cell(2:npts+1,1) = num2cell(Solution.t/3600);
wrt_cell(2:npts+1,2:nvar+1) = num2cell(Solution.y');
status = xlswrite('sim_results',wrt_cell,sheet);
cd(curr_dir);
return