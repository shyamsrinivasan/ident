%% Get Data & Obtain Initial Solution
%Build 1.1
%New mass balance for gene expression <=> No more Elasticity matrix
%Separated basal expression from regulated expression
%No Self Regulation
%Does not account for input metabolites 

%April 17 Accounts for input metabolites (most of it!)
%Does not work: galS, galR, mlc, arcA

%April 30 Works for all genetic regulons using newloglinear and linlogfunc1
%Need to add metabolites
%Solve initially -> Use SS values as initial values to test any subsequent
%changes to parameters/system

%July 27 Modified Input Metabolite Function and streamlined code
%Initialization of parameters as vectors as opposed to scalars
% => Parameters are independent for each gene
%Probability of using a different model for high-dimensional interactions
%%

clc
clear all
%load('D:\My Documents\Courses\CHE 1125H\Project\Ecoli Global1.mat');
%load('D:\My Documents\Courses\CHE 1125H\Project\Reduced Ecoli Model.mat');
load('C:\Users\shyam\Documents\Courses\CHE 1125H\Project\TRN Model\Reduced Ecoli Model.mat');
%load('D:\My Documents\Courses\CHE 1125H\Project\RNDNUM.mat');

%Assigning Input Metabolites
inputmetab = {'GLCxt';'O2xt'};
MetabInConc = [10;10];
[MetabIn] = changeinput(model,inputmetab,MetabInConc);

%Generate Pseudo Genes
%PGenes = model.PGMat*MetabIn;
%PGenes(PGenes>0)=1;
PGenes = [];

%Generate Vectors
[TF,DNA,G,PGenes,EGene,EPSGene]=GenerateVectors(model,PGenes);

%Parameter Initialization  # of parameters: 8 

ngenes = size(G,1);
ntf = size(TF,1);
parameters = struct([]);

parameters(1).basaltr = zeros(ngenes,1);
parameters.regtr = zeros(ngenes,1);
parameters.hillcoeff = zeros(ngenes,1);
parameters.actcoeff = zeros(ngenes,1);
parameters.repcoeff = zeros(ngenes,1);



%basal transcription rate, kb
parameters.basaltr(parameters.basaltr==0) = 0.0133;         
%stimulated transcription rate, ks
parameters.regtr(parameters.regtr==0) = 0.0016;           
%Hill coefficient n for repression
parameters.hillcoeff(parameters.hillcoeff==0) = 2;
%beta
parameters.actcoeff(parameters.actcoeff==0) = 0.008;        
%theta
parameters.repcoeff(parameters.repcoeff==0) = 0.08;        
%Translation rate
parameters.ktranslation = 0.0006;
%Gene Decay
parameters.genedecay = 0.00016;
%protein decay        
parameters.protdecay = 0.00016;

parameters.beta1 = parameters.actcoeff;           %Metabolite Activation
parameters.theta1 = parameters.repcoeff;          %Metabolite Repression

%Solving to Establish Initial SS   

file = 'C:\Users\shyam\Documents\Courses\CHE 1125H\Project\Results.xlsx';
T = {'indigo';'magenta';'cyan';'red';'green';'blue';'Lime';'Aqua';'Crimson';'Black'};

%[time,dG,ngenes,ntf] = testlinear(parameters,TF,DNA,G,PGenes,EGene,EPSGene,MetabIn);
%G = dG(end,1:ngenes)';
%TF = dG(end,ngenes+1:end)';

%writetofile(file,G,'Sheet1');
%writetofile(file,TF,'Sheet2');

%[time,dG] = testlinear(P,TF,DNA,Gene1,PGenes,EGene,EPSGene,MetabIn);
%[time,dG] = newloglinear(P,TF,DNA,Gene1,PGenes,EGene,EPSGene,MetabIn);

clc
%Parameter Adjustment for individual Genes for local sensitivity analysis
val = [8E-1;8E-2;8E-3;8E-4;8E-5;8E-6];
sheet = {'Sheet3','Sheet4','Sheet5','Sheet6','Sheet7','Sheet8'};
nval = length(val);
%inputgene = {'fdnG','pflA'};
%inputgene = {'fdnG';'focA';'fumB';'pflA';'sdhA'};
inputgene = {'aceA';'fadR';'iclR';'mlc';'ptsG'};
ngeneinput = length(inputgene);
T = zeros(ngenes+ntf,1);
igene = 1;
while igene <= ngeneinput
    geneindx = strcmp(inputgene{igene},model.GeneName);
    T(geneindx) = 1;
    igene = igene + 1;
end
T = logical(T');
Solution.Gss = cell(nval,1);
Solution.Gpert = cell(nval,1);
Solution.foldchange = cell(nval,1);

flag = 0;
for ival = 1:nval
    parname = {'repcoeff'};
    %parname = {'actcoeff'};
    parval = {val(ival);val(ival);val(ival);val(ival);val(ival)};
    %[parameters,flag] = changeparameter(model,parameters,gene,parname,parval)
    [parameters] = changeparameter(model,parameters,inputgene,parname,parval);
    
    if ~flag
        inputmetab = {'GLCxt';'O2xt'};
        MetabInConc = [10;10];
        [MetabIn] = changeinput(model,inputmetab,MetabInConc);

        [time,dG,ngenes,ntf,model] = testlinear(parameters,TF,DNA,G,PGenes,EGene,EPSGene,MetabIn,model);
        Solution.Gss{ival} = dG(end,T)';
        Solution.G = dG(end,1:ngenes)';
        %writetofile(file,G,'Sheet1');
        G = Solution.G;
    
        inputmetab = {'GLCxt';'Acxt'};
        MetabInConc = [0;10];
        [MetabIn] = changeinput(model,inputmetab,MetabInConc,MetabIn);
        [time,dG,~,~,model] = testlinear(parameters,TF,DNA,G,PGenes,EGene,EPSGene,MetabIn,model);
        Solution.Gpert{ival} = dG(end,T)';
        Solution.newG = dG(end,1:ngenes)';
    
        Solution.allFC = log2(Solution.G./Solution.newG);
        %newGE = dG(end,T);    
    
    
        Solution.foldchange{ival} = log2(Solution.Gss{ival}./Solution.Gpert{ival});
        Solution.InputGene = model.GeneName(T);
        
    end
    
    
  
end

%%

for ival = 1:nval
    
    
    inputmetab = {'O2xt'};
    MetabInConc = [0];
    [MetabIn] = changeinput(model,inputmetab,MetabInConc,MetabIn);
    [time,dG] = testlinear(parameters,TF,DNA,G,PGenes,EGene,EPSGene,MetabIn);
    GE = dG(end,1:ngenes)';
    
    
    ngeneinput = length(inputgene);
    T = zeros(ngenes+ntf,1);
    igene = 1;
    while igene <= ngeneinput
        geneindx = strcmp(inputgene{igene},model.GeneName);
        T(geneindx) = 1;
        igene = igene + 1;
    end
    newGE = dG(end,logical(T'));    
    writetofile(file,newGE,sheet{ival});

end







%%

i=10;
%while i<=10
   
    inputmetab = {'O2xt'};
    MetabInConc = [0];
    [MetabIn] = changeinput(model,inputmetab,MetabInConc,MetabIn);
    %[time,dG] = newloglinear(P,TF,DNA,Gene1,PGenes,EGene,EPSGene,MetabIn);
    [time,dG] = testlinear(parameters,TF,DNA,G,PGenes,EGene,EPSGene,MetabIn);
    GE = dG(end,1:ngenes)';
    %PE = dG(end,ng+1:end)';
    writetofile(file,GE,'Sheet2');
    %writetofile(file,TF,'Sheet4');
   
    
    %Fold Change
    %G_FC = foldchange(G,GE);
    %P_FC = foldchange(TF,PE);
    %writetofile(file,G_FC,'Sheet5');
    %writetofile(file,P_FC,'Sheet6');
    
    time = time/3600;
    %dG = log(dG);
    %subplot(2,4,1);
    index = strncmp('galS',model.GeneName,4);
    subplot(2,4,1);
    index=find(index)+80;
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('galS Expression');
    xlabel('Time (h)')
    hold on
    
    index = strncmp('galE',model.GeneName,4);
    subplot(2,4,2);
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('galE Expression');
    xlabel('Time (h)')
    hold on
    
    index = strncmp('galR',model.GeneName,4);
    subplot(2,4,4);
    index=find(index)+80;
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('galR Expression');
    xlabel('Time (h)')
    hold on
    
    
    
    index = strncmp('arcA',model.GeneName,4);
    subplot(2,4,3);
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('arcA Expression');
    xlabel('Time (h)')
    hold on 
    
    
 
    index = strncmp('mdh',model.GeneName,3);
    subplot(2,4,5);
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('mdh Expression');
    xlabel('Time (h)')
    hold on
    
    index = strncmp('sucA',model.GeneName,4);
    subplot(2,4,8);
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('sucA Expression');
    xlabel('Time (h)')
    hold on

    index = strncmp('fdnG',model.GeneName,4);
    subplot(2,4,6);
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('fdnG Expression ');
    xlabel('Time (h)')
    hold on
    
    index = strncmp('pdhR',model.GeneName,4);
    subplot(2,4,7);
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('pdhR Expression');
    xlabel('Time (h)')
    hold on
        
    i=i+1;
%end
%ylabel('ln(mRNA)');
legend;
    %% Plot Calls for Various Genes
    subplot(2,4,4);
    index = strncmp('fadR',model.GeneName,4);
    index=find(index)+139;
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('fadR Expression');
    xlabel('Time (h)')
    hold on
    
    index = strncmp('iclR',model.GeneName,4);
    index=find(index)+139;
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('iclR Expression');
    xlabel('Time (h)')
    hold on
    
    index = strncmp(' aceA',model.GeneName,5);
    subplot(2,4,2);
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('aceA Expression');
    xlabel('Time (h)')
    hold on
    
   index = strncmp(' adhE',model.GeneName,5);
    subplot(2,4,3);
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('adhE Expression');
    xlabel('Time (h)')
    hold on
    
    
    
    index = strncmp(' aceB',model.MetabNames,5);
    subplot(5,1,4);
    plot(time,PGenes(index));
    ylabel('aceB Expression');
    
    subplot(nofplots,1,2);
    index = strncmp(' sdhA',model.GeneName,5);
    plot(time,dGeneinit(:,index));
    xlabel('sdhA Expression');
    
    subplot(nofplots,1,3);
    index = strncmp('pflA',model.GeneName,4);
    plot(time,dGeneinit(:,index));
    xlabel('pflA Expression');
    
    subplot(2,4,4);
    index = strncmp(' fnr',model.GeneName,4);
    index=find(index)+139;
    plot(time,dG(:,index),'color',rgb(T(i)));
    ylabel('fnr Expression');
    xlabel('Time (h)')
    hold on
    
    %subplot(nofplots,1,4);
    %index = strncmp(' fumB',model.GeneName,5);
    %plot(time,dGeneinit(:,index));
    %xlabel('fumB Expression');
 %UnregulTED Genes
    %K = {5,10,13,14,22,27,28,30,31,35,36,45,53,54,55,63,64,65,66,67,68,73,75,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,100,102,103,106,107,108,109,111,115,123,124,131,134,135,136,137,138,139};
%%



%%
clear all
[num,col_heads] = xlsread('D:\My Documents\Courses\CHE 1125H\Project\RMatrix.xlsx','Sheet3');
FCM_Oxg = sign(num(:,1));
FCM_Act = sign(num(:,2));
FCM_Gly = sign(num(:,3));
FCE_Oxg = sign(num(:,4));
FCE_Act = sign(num(:,5));
FCE_Gly = sign(num(:,6));

M1 = zeros(80,3);
M1(FCM_Oxg==FCE_Oxg,1)=1;
M1(FCM_Act==FCE_Act,2)=1;
M1(FCM_Gly==FCE_Gly,3)=1;
M1 = [M1;sum(M1,1);sum(M1,1)*100/80]; 



%%
%Miscellaneous Commands
      
%[TF,DNA1,DNA2,Gene1,Gene2] = GenerateRegulatedTF(model,MetabIn);
%[ETF, EDNA1,EDNA2,EGene,TF,DNA1,DNA2,Gene1,Gene2] =
%GenerateE(E,TF,DNA1,DNA2,Gene1,Gene2);

%[ETF, EDNA,EGene,TF,DNA,Gene] = GenerateE(E,TF,DNA,dGenefinal);
%[time,dGene1] = loglinear(e,ETF,EDNA,EGene,TF,DNA,Gene);
    

%E = xlsread('D:\My Documents\Courses\CHE 1125H\Project\RMatrix.xlsx','EPSGene');

%[num,model.TFname] = xlsread('C:\Users\B4-x32\Documents\Courses\CHE 1125H\Project\RMatrix.xlsx','

%xlswrite('filename',matrix,'sheet','rangelb:rangeub');
%xlswrite('D:\My Documents\Courses\CHE 1125H\Project\Data.xlsx',EGene,'Sheet1');
%t = cell2mat(sol.e);
%t = num2str(t);

%sol.dG{i,1} = real(dG);
%sol.t{i,1} = time;
%sol.e{i,1} = e;


%%
%Parameter Estimation & Metabolic Network
%%Random Number Generation
%s_qmc=sobolset(1,'Skip',100);
%z_qmc=net(s_qmc,nofgenes);

%Sampling Algorithms
%for i=1:N
%random number = 
%Proposal dustribution = q() to determine the next point in the markov
%chain
%q() should be as close as possible to p()(the original probability density
%function that you are trying o sample from

%From testlinear function

switch flag(i)
              case 0%No R
                  dV(i) = KB(i)*(DNA(i))-kd*V(i);
          
              case 1%ER>0
                  dV(i) = KS(i)*sum((TF(EGene(i,:)>0)*beta(i))./...
                          (beta(i)+TF(EGene(i,:)>0)))-KD*V(i);
                  
              case 2%ER<0
                  dV(i) = KS(i)*sum(1/...
                          (1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))-KD*V(i);
                  
              case 3%ER>0 && ER<0
                  dV(i) = KS(i)*(sum(1/(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))+...
                          sum((TF(EGene(i,:)>0)*beta(i))./(beta(i)+TF(EGene(i,:)>0))))...
                          -KD*V(i);
                  
              case 4%MR>0
                  dV(i) = KS(i)*sum((M(EPSGene(i,:)>0)*beta1(i))./(beta1(i)+M(EPSGene(i,:)>0)))-...
                          KD*V(i);
                  
              case 5%MR<0
                  dV(i) = KS(i)*(sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i))))-...
                          KD*V(i);%M(EPSGene(i,:)<0)
                  
              case 6%MR>0 && MR<0
                  dV(i) = KS(i)*(sum((M(EPSGene(i,:)>0)*beta1(i))./(beta1(i)+M(EPSGene(i,:)>0)))+...
                          sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i))))-...
                          KD*V(i);%M(EPSGene(i,:)<0).
                  
              case 7%ER>0 && MR>0
                  dV(i) = KS(i)*(sum((TF(EGene(i,:)>0)*beta(i))./(beta(i)+TF(EGene(i,:)>0))))*...
                          (sum((M(EPSGene(i,:)>0)*beta1(i))./(beta1(i)+M(EPSGene(i,:)>0))))-...
                          KD*V(i); 
                  
              case 8%ER>0 && MR<0
                  dV(i) = KS(i)*(sum((TF(EGene(i,:)>0)*beta(i))./(beta(i)+TF(EGene(i,:)>0))))*...
                          (sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i))))-...
                          KD*V(i);%M(EPSGene(i,:)<0).
                  
              case 9%ER<0 && MR>0
                  dV(i) = KS(i)*sum(1./(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))*...
                          sum((M(EPSGene(i,:)>0)*beta1(i))./(beta1(i)+M(EPSGene(i,:)>0)))-...
                          KD*V(i);
                  
              case 10%ER<0 && MR<0
                  dV(i) = KS(i)*sum(1/(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))*...
                          sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i)))-...
                          KD*V(i);%TF(EGene(i,:)<0). M(EPSGene(i,:)<0).
                  
              case 11%ER>0 && ER<0 && MR<0
                  dV(i) = KS(i)*(sum(1./(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))+...
                          sum((TF(EGene(i,:)>0)*beta(i))./(beta(i)+TF(EGene(i,:)>0))))*...
                          sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i)))-...
                          KD*V(i);%M(EPSGene(i,:)<0).
                  
              case 12%ER>0 && ER<0 && MR>0
                  dV(i) = KS(i)*(sum(1./(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))+...
                          sum((TF(EGene(i,:)>0)*beta(i))./(beta(i)+TF(EGene(i,:)>0))))*...
                          sum((M(EPSGene(i,:)>0)*beta1(i))./(beta1(i)+M(EPSGene(i,:)>0)))-...
                          KD*V(i);
                  
              case 13%ER>0 && MR<0 && MR>0
                  dV(i) = KS(i)*(sum((M(EPSGene(i,:)>0)*beta1(i))./(beta1(i)+M(EPSGene(i,:)>0)))+...
                          sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i))))*...
                          sum((TF(EGene(i,:)>0)*beta(i))./(beta(i)+TF(EGene(i,:)>0)))-...
                          KD*V(i);%M(EPSGene(i,:)<0).
                  
              case 14%ER<0 && MR<0 && MR>0
                  dV(i) = KS(i)*(sum((M(EPSGene(i,:)>0)*beta1(i))./(beta1(i)+M(EPSGene(i,:)>0)))+...
                          sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i))))*...
                          sum(1./(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))-...
                          KD*V(i);%M(EPSGene(i,:)<0).
          end


                            
                        
                   
    
    
    