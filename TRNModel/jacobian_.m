function J = jacobian_(ng,eq_pt,data,model,FBAmodel)
nreg = sum(ng(2:end));
mRNA = eq_pt(1:ng(1));
regulator = eq_pt(ng(1)+1:ng(1)+nreg);
Mxt = eq_pt(ng(1)+ng(2)+ng(3)+1:sum(ng));
alpha = 0.01;%transcription rate
pdecay = data.pdecay;
mdecay = data.mdecay;
trate = 1;
b = 20;
Rmax = 1e-5;
Km = 50;
nvar = sum(ng)+1;
intS = FBAmodel.S(1:ng(3),:);
uptakeS = FBAmodel.M*FBAmodel.S(1:ng(3),FBAmodel.Vic_exind);
SSval = model.SSval;

% V_ind = [FBAmodel.Vin_ind,FBAmodel.Vout_ind];
% uptakeS = FBAmodel.M*FBAmodel.S(1:ng(3),V_ind);


J = zeros(nvar,nvar);
% J(1:ng(1),1:ng(1)) = sparse(1:ng(1),1:ng(1),-mdecay,ng(1),ng(1));




% 
% J(end,ng(1)+ng(2)+ng(3)+1:sum(ng)) = Rmax./(Km+Mxt)-Rmax.*Mxt./(Km+Mxt);
% J(1:ng(1),end) = (1+b)/b;

metab_full = SSval(ng(1)+ng(2)+1:sum(ng)).*(1+regulator(ng(2)+1:sum(ng)-ng(1)));

%Indices
%1:ng(1) - mRNA
%ng(1)+1:ng(1)+ng(2) - Proteins
%ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3) - Internal Metabolites
%ng(1)+ng(2)+ng(3)+1:sum(ng) - External Metabolites
%end - RNAP

J = zeros(nvar,nvar);
%mRNA
J(1:ng(1),1:ng(1)) = sparse(1:ng(1),1:ng(1),-mdecay,ng(1),ng(1));
J(1:ng(1),end) = sparse(1:ng(1),1,alpha*SSval(end)./SSval(1:ng(1)),ng(1),1);
%Protein
JPG = trate*SSval(1:ng(1))./SSval(ng(1)+1:ng(1)+ng(2));
JPG = diag(JPG);
JPP = eye(ng(2),ng(2));
JPP(JPP==1) = -pdecay;
JP = [JPG JPP];
J(ng(1)+1:ng(1)+ng(2),1:ng(1)+ng(2)) = JP;
%Intracellular Metabolite
[m_ind,reg_ind] = find(intS);
l_ind = sub2ind(size(intS),m_ind,reg_ind);
Jl_ind = sub2ind(size(J),m_ind+ng(1)+ng(2),reg_ind+ng(1));
JiMP1 = intS(l_ind).*SSval(reg_ind+ng(1))./SSval(m_ind+ng(1)+ng(2));
J(Jl_ind) = 0.1.*JiMP1;
JiMiM = eye(ng(3));
JiMiM(JiMiM==1) = -0.2*SSval(end)*(1+eq_pt(end));
J(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3),ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3)) = JiMiM;
%Extracellular Metabolite
[mx_ind,regx_ind] = find(uptakeS);
lx_ind = sub2ind(size(uptakeS),mx_ind,regx_ind);
Jlx_ind = sub2ind(size(J),mx_ind+ng(1)+ng(2)+ng(3),regx_ind+FBAmodel.n_rxn+ng(1));
JxMP1 = -uptakeS(lx_ind).*SSval(regx_ind+FBAmodel.n_rxn+ng(1))./SSval(mx_ind+ng(1)+ng(2)+ng(3));
J(Jlx_ind) = JxMP1;
JxMxM = eye(ng(4));
JxMxM(JxMxM==1) = -1;
J(ng(1)+ng(2)+ng(3)+1:sum(ng),ng(1)+ng(2)+ng(3)+1:sum(ng)) = JxMxM;
%rnap
Jrnap1 = ([1e-5;0].*SSval(ng(1)+ng(2)+ng(3)+1:sum(ng)))./([50;1]+metab_full(ng(3)+1:ng(3)+ng(4)));
Jrnap2 = (metab_full(ng(3)+1:ng(3)+ng(4)).*[1e-5;0]).*(SSval(ng(1)+ng(2)+ng(3)+1:sum(ng)));
JRM = (1/SSval(end))*(Jrnap1-Jrnap2);
J(end,ng(1)+ng(2)+ng(3)+1:sum(ng)) = JRM;
J(end,end) = -pdecay;




% J(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3),ng(1)+1:ng(1)+ng(2)) = FBAmodel.S(1:ng(3),:)*model.maxReg(1:ng(2));
% J(ng(1)+ng(2)+ng(3)+1:sum(ng),ng(1)+1:ng(1)+ng(2)) = -FBAmodel.M*FBAmodel.S(1:ng(3),V_ind)*model.maxReg(V_ind);
 

return







