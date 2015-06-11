function [time,dGene] = loglinear(e,TF,DNA1,DNA2,Gene1,Gene2,Gene3,ETF,EDNA1,EDNA2,EGene,EPSGene)
%Version 1.0a
% Using log G instead of G in the decay terms for Genes and TFs. The
% profiles are strakly different from that witnessed when G was used
      
    regGenes = size(Gene1,1);
    tfGenes = size(Gene2,1);
   % psgenes = size(Gene3,1);
    tfs = size(TF,1);
    
    V = [Gene1;Gene2;TF];
    kdecay = 0.16;
    ktranslation = 0.01;
    Vinit = V';
    freq_const = 0.89;
    
    
    [time, dG] = ode23(@linlogfunc,[0 200],Vinit);
    dGene = dG;
    %time = T;
    %dGene = vratio;
    %fig(T,vratio,4);
    
    
    function dV = linlogfunc(t,V)
        dV = ones(regGenes+tfGenes+tfs,1);
        G = V(1:regGenes);
        TFG = V(regGenes+1:regGenes+tfGenes);
        TF = V(regGenes+tfGenes+1:regGenes+tfGenes+tfs);
       
%Original Equation

        dV(1:regGenes) = e*(1+ETF*log(TF)+EDNA1*log(DNA1)+EGene*log(G)+EPSGene*log(Gene3))-kdecay*G;
        dV(regGenes+1:regGenes+tfGenes) = e*(1+EDNA2*log(DNA2))-kdecay*TFG;
        dV(regGenes+tfGenes+1:regGenes+tfGenes+tfs) = ktranslation*TFG - kdecay*TF;
        
        
        %for i=1:size(EGene,1)
         %   for j=1:size(EGene,2)
         %       if TF(j)>0 && EGene(i,j)==-1
         %           dY(i)=0;
         %       else 
         %           dY(i) = e*(1+ETF(i,:)*log(TF)+EDNA(i,:)*log(DNA)+EGene(i,:)*log(Gene));
         %       end
         %   end
        %end
    end


    function fig(time,vratio,plotno)
        %nofgenes = size(dGene,1);
        %timestep = size(time,1);
        for k=1:plotno
            subplot(plotno,1,k)
            plot(time,vratio(:,k));
        end
       
  
    end
        
end