function [kcat_ratio,kcat_fwd,kcat_bkw] =...
         haldane_kcat(model,KVl,subs,pruds,irxn,pvector)
% 
% low = 0.1;
% high = 2000;
% nsamp = 1000;

Km = pvector.K;
Keq = model.Keq;
Kcat = model.Kcat;
Vss = model.Vss;
     
s_sub = model.S(subs,irxn);
s_prd = model.S(pruds,irxn);

lnkcat = log(Keq(irxn)) -...
         sum(s_sub.*log(Km(subs,irxn))) + sum(s_prd.*log(Km(pruds,irxn)));
kcat_fwd = model.Kcat(irxn);
kcat_bkw = exp(log(Kcat(irxn)) - lnkcat/2);  

if Vss(irxn) < 0
    if kcat_bkw < kcat_fwd
        kcat_lb = exp(log(KVl(irxn)) + lnkcat/2);
        kcat_ub = 2*exp(log(KVl(irxn)) + lnkcat/2);        
        kcat_bkw = kcat_lb + (kcat_ub-kcat_lb).*rand(1,1);
%         kcat_fwd = exp(log(KVl(irxn)) + lnkcat/2);
    end
elseif Vss(irxn) > 0
    if kcat_bkw > kcat_fwd
        kcat_lb = 0;
        kcat_ub = exp(log(KVl(irxn)) + lnkcat/2);        
        kcat_bkw = kcat_lb + (kcat_ub-kcat_lb).*rand(1,1);
%         kcat_fwd = exp(log(KVl(irxn)) + lnkcat/2);
    end
else
    kcat_fwd = 0;
    kcat_bkw = 0;
end
kcat_ratio = exp(lnkcat);%kcat_fwd/kcat_bkw


% if ~isnan(KVl(irxn))
%     kcat_fwd = exp(log(KVl(irxn)) + lnkcat/2);
%     kcat_bkw = exp(log(KVl(irxn)) - lnkcat/2);
% else%sample KVl
%     if model.rev(irxn) && model.Vss(irxn) < 0
%         KVl = low + (high-low).*rand(nsamp,1);
%         kcat_fwd = exp(log(KVl(irxn)) + lnkcat/2);
%         kcat_bkw = exp(log(KVl(irxn)) - lnkcat/2);
%         pvector.kcat_ratio = kcat_ratio;
%         pvector.kcat_fwd = kcat_fwd;
% %         flag = checkVmax(pvector);
% %         while kcat_bkw < kcat_fwd && flag < 0
% %             KVl = low + (high-low).*rand(nsamp,1);
% %             kcat_fwd = exp(log(KVl(irxn)) + lnkcat/2);
% %             kcat_bkw = exp(log(KVl(irxn)) - lnkcat/2);
% %             pvector.kcat_ratio = kcat_ratio;
% %             pvector.kcat_fwd = kcat_fwd;
% %             flag = checkVmax(pvector);            
% %         end       
%         
%     end
% end



