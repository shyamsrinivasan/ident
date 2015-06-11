function [varargout] =...
globalSens_serial(model,FBAmodel,ng,defparval,varname,saveData,varargin)
Amat = varargin{1};
Bmat = varargin{2};
Cmat = varargin{3};
varargout = cell(3,1);
yA = zeros(sum(ng)+1,nsampl);
yB = zeros(sum(ng)+1,nsampl);
yC = zeros(sum(ng)+1,nsampl,Ns);
for jsampl = 1:nsampl
    fprintf('Sample: %d\n',jsampl);
    %Assign parameters from A
    fprintf('Calculating from Matrix A\n');
    model = sensAssign(model,Amat,jsampl);
    %Simulate Model
    [model,batch,solverP] = initializeModel(model,36000); 
    %Solve Model ODE 
    initSolution = struct([]);
    [~,~,~,finalSS] =...
    runSolver(model,FBAmodel,batch,defparval,ng,...
                solverP,initSolution,...
                varname,saveData);  
    yA(:,jsampl) = finalSS.y;
    close all         
    fprintf('Calculating from Matrix B\n');
    %Assign parameters from B
    model = sensAssign(model,Bmat,jsampl);
    %Simulate Model
    [model,batch,solverP] = initializeModel(model,36000); 
    %Solve Model ODE 
    initSolution = struct([]);
    [~,~,~,finalSS] =...
    runSolver(model,FBAmodel,batch,defparval,ng,...
                solverP,initSolution,...
                varname,saveData);  
    yB(:,jsampl) = finalSS.y;
    close all
    fprintf('Calculating from Matrix C\n');
    for j = 1:Ns
        %Assign parameters from C
        model = sensAssign(model,Cmat(:,:,j),jsampl);
        %Simulate Model
        [model,batch,solverP] = initializeModel(model,36000); 
        %Solve Model ODE 
        initSolution = struct([]);
        [~,~,~,finalSS] =...
        runSolver(model,FBAmodel,batch,defparval,ng,...
                    solverP,initSolution,...
                    varname,saveData);  
        yC(:,jsampl,j) = finalSS.y;
        close all 
    end            
end
varargout{1} = yA;
varargout{2} = yB;
varargout{3} = yC;

return