function [time,dGene] = newloglinear(P,TF,DNA,Gene1,PGenes,EGene,EPSGene,M)
%linlogfunc1 - Highest accuracy/approximation closest to reality
%Changes to this function will be tested in testlinear()
      
    ngenes = size(Gene1,1);
    ntf = size(TF,1);

    KD = zeros(ngenes,1);    %9 Parameters 
    KD(KD==0)=P(1);          %Gene Regulated Decay
    kd = P(2);               %Gene Basal Decay
    KB = zeros(ngenes,1);
    KB(KB==0)=P(3);          %Basal Expression
    KS = zeros(ngenes,1);
    KS(KS==0)=P(4);          %Regulated Expression
    ktranslation = P(5);     %Protein Translation
    kdecay = P(6);           %Protein Decay
    beta = P(7);             %beta -> Activation
    theta = P(8);            %theta -> Repression
    n = P(9);                %Repression Exponent
        
    V = [Gene1;TF];
    Vinit = V';
    [time, dG] = ode23(@linlog,[0 36000],Vinit);
    dGene = dG;
    function dV = linlog(t,V)
       dV = ones(ngenes+ntf,1);
       G = V(1:ngenes);
       TF = V(ngenes+1:end);
       
       for i=1:ngenes
          Ractivated = 0;
          Rinhibited=0;
          RMactivated=0;
          RMinhibited=0;
          flag=0;
          if ~isempty(find(EGene(i,:)~=0,1)) 
             if ~isempty(find(EGene(i,:)<0,1)) && ~isempty(find(EGene(i,:)>0,1))
                dV(i) = KS(i)*(sum(TF(EGene(i,:)<0)./(1+(TF(EGene(i,:)<0)/theta).^n))+...
                        sum((TF(EGene(i,:)>0)*beta)./(beta+TF(EGene(i,:)>0))))...
                        -KD(i)*V(i);
             else if ~isempty(find(EGene(i,:)<0,1))
                     dV(i) = KS(i)*sum(TF(EGene(i,:)<0)./...
                             (1+(TF(EGene(i,:)<0)/theta).^n))-KD(i)*V(i);
                     flag=1;
                 else if ~isempty(find(EGene(i,:)>0,1))
                         dV(i) = KS(i)*sum((TF(EGene(i,:)>0)*beta)./...
                                 (beta+TF(EGene(i,:)>0)))-KD(i)*V(i);
                         flag=1;
                     end
                 end
             end
          else
              dV(i) = KB(i)*(DNA(i))-kd*V(i);
          end                 
       end
       dV(ngenes+1:end) = ktranslation*V(1:ngenes)-kdecay*V(ngenes+1:end);
    end
 
end        
