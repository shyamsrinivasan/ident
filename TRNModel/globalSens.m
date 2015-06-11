% % globalSens(model,nsampl,parName,parScale)
clc
% cd('C:\Users\shyam\Documents\Courses\CHE 1125 Project'); 
load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Kinetic Model\N2m.mat');
% load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\integrated_space.mat');
addpath(genpath('C:\Users\shyam\Documents\MATLAB\sundialsTB'));
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model'));
addpath(genpath('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Other Matlab Central Files'));
trnfname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\N2trn.txt';
parameter = ensb.model1;
[trnmodel,FBAmodel,defparval] = Tmodel(trnfname,FBAmodel,parameter,variable);
ng = [trnmodel.nt_gene;...
      trnmodel.nt_prot+trnmodel.nr_cmplx-length(trnmodel.PMind_R)-1;...      
      length(trnmodel.PMind_R);...
      trnmodel.nt_metab-length(trnmodel.PMind_R)-trnmodel.next_metab-length(trnmodel.bm_ind);...
      0;...%length(trnmodel.RegCMPLX);...
      trnmodel.next_metab]; 
%Sensitivity of model to parameters
model = trnmodel;
nsampl = 10;
parScale = [1e-4 20;...%Coefficient
            500 5000;...%brate
            0.01 1;...%ptrate
            0 0;...%
            0 0;...
            0.1 1.0;...%gmax
            0.005 0.5;...%Vuptake
            1   30000];%Kcat
parList = {'Coefficient','brate','ptrate','Kb','Kub','gmax','Vuptake'};
[model,Amat,Bmat,Cmat,Ns] = globalSensMatrix(model,nsampl,parScale);
%Assign parameters from matrices for model computation
%Each column of each matrix is new set of parameters
varname = {'A';'P';'P4';'P1'}; 
saveData.filename = '';
saveData.dirname =...
'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Results\TRN Model version 7';
%Parallel Implementaion of global Snesitivity Analysis
cluster = parcluster;
j = createJob(cluster);
% inputs = {model,FBAmodel,ng,Amat};
createTask(j,@globalSens_parallel,1,{{model,FBAmodel,ng,defparval,varname,saveData,Amat},...
                                     {model,FBAmodel,ng,defparval,varname,saveData,Bmat}});
for jcol = 1:Ns
    C = Cmat(:,:,jcol);    
    createTask(j,@globalSens_parallel,1,{model,FBAmodel,ng,defparval,varname,saveData,C});
%     submit(j);
end
submit(j);
wait(j);
data = fetchOutputs(j);
delete(j);
yA = data{1};
yB = data{2};
yC = zeros(sum(ng)+1,nsampl,Ns);
for jcol = 1:Ns
    yC(:,:,jcol) = data{2+jcol};
end
[S,ST] = SensitivityIndex(yA,yB,yC,Ns);
% S = zeros(sum(ng)+1,Ns);
% for j = 1:Ns
%     f02 = sum(yA,2);
%     nr = diag(yA*yC(:,:,j)') - f02;
%     dr = diag(yA*yA') - f02;
%     S(:,j) = nr./dr;
% end
% j1 = createCommunicatingJob(cluster,'Type','pool');
% j2= createCommunicatingJob(cluster,'Type','pool');
% j3 = createCommunicatingJob(cluster,'Type','pool');
% j.NumWorkersRange = [1 2];
%wait(j);
% data = fetchOutputs(job);
% delete(job);

% load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model\sensData.mat');






