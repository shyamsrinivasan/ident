%Parallel Implementaion of global Snesitivity Analysis
function [y] = globalSens_parallel(trnmodel,FBAmodel,ng,defparval,varname,saveData,Mat)
nsampl = size(Mat,1);    
y = zeros(sum(ng)+1,nsampl);
for jsampl = 1:nsampl    
    fprintf('Sample: %d\n',jsampl);
    %Assign parameters from Mat
%     fprintf('Calculating from Matrix A\n');
    M = Mat(jsampl,:);
    model = sensAssign(trnmodel,M);
    %Simulate Model
    [model,batch,solverP] = initializeModel(model,36000); 
    %Solve Model ODE 
    initSolution = struct([]);
    [~,~,~,finalSS] =...
    runSolver(model,FBAmodel,batch,defparval,ng,...
                solverP,initSolution,...
                varname,saveData);  
    y(:,jsampl) = finalSS.y;
    close all         
end
end