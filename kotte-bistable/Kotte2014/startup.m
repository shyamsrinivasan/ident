% startup.m
if ~exist('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel')
    status = 2;
    fprintf('\nLinux System\n');
else 
    status = 1;
    fprintf('\nWindows System\n');
end

if status == 1
    addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
    addpath('C:\Users\shyam\Documents\MATLAB\CCFM_manifolds\CCFM_manifolds\functions\');    
elseif status == 2
    addpath('/home/shyam/Documents/MATLAB/Code');
    addpath(genpath('/home/shyam/Documents/MATLAB/Code/KineticModel'));
    addpath(genpath('/home/shyam/Documents/MATLAB/Code/CBM'));
    addpath(genpath('/home/shyam/Documents/MATLAB/Code/plots/'));   
    addpath('/home/shyam/Documents/MATLAB/Toolbox/CCFM_manifolds/CCFM_manifolds/functions');    
    addpath(genpath('/home/shyam/Documents/MATLAB/Toolbox/MatlabCentral'));
    addpath(genpath('/home/shyam/Documents/MATLAB/Toolbox/zz_ADMAT-2.0'));
    addpath('/home/shyam/Documents/MATLAB/mat_files');
    addpath('/home/shyam/Documents/MATLAB/Toolbox');
    cd('/home/shyam/Documents/MATLAB/Toolbox/Cl_matcont5p4/Cl_matcont5p4');
    init
    cd('/home/shyam/Documents/MATLAB/Code/KineticModel/Kotte2014');
    addpath(genpath('/home/shyam/Documents/MATLAB/Toolbox/Cl_matcont5p4'));
end
    
% try
%     
% catch
%     
%     
%     addpath('/home/shyam/Documents/MATLAB/Toolbox/CCFM_manifolds/CCFM_manifolds/functions');
%     addpath(genpath('/home/shyam/Documents/MATLAB/Toolbox/Cl_matcont5p4'));
%     addpath(genpath('/home/shyam/Documents/MATLAB/Toolbox/MatlabCentral'));
%     addpath(genpath('/home/shyam/Documents/MATLAB/Toolbox/zz_ADMAT-2.0'));
%     addpath('/home/shyam/Documents/MATLAB/mat_files');
%     addpath('/home/shyam/Documents/MATLAB/Toolbox');
%     rxfname = '/home/shyam/Documents/MATLAB/Code/KineticModel/Kotte2014/Kotte2014.txt';
%     cnfname = '/home/shyam/Documents/MATLAB/Code/KineticModel/Kotte2014/Kotte2014C.txt';
%     fprintf('\nLinux System\n');
% end