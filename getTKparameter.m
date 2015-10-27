function pvec = getTKparameter(model,pvec,mc,Vex)


h2o = find(strcmpi(model.mets,'h2o[c]'));
pic = find(strcmpi(model.mets,'pi[c]'));
pie = find(strcmpi(model.mets,'pi[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
he = find(strcmpi(model.mets,'h[e]'));

for irxn = 1:length(Vex)
    if model.Vss(Vex(irxn))~=0
        %kinetics - substrate and product
        sbid = logical(model.S(:,Vex(irxn))<0);
        prid = logical(model.S(:,Vex(irxn))>0);
        if any(sbid) && any(prid)
            if ~strcmpi(model.rxns{Vex(irxn)},'h2ot')
                sbid(h2o) = 0;
            end
            if ~strcmpi(model.rxns{Vex(irxn)},'PIt2r')
                sbid([pic pie hc he]) = 0;
                prid([pic pie hc he]) = 0;
            else
                sbid([hc he]) = 0;
                prid([hc he]) = 0;
            end
            
            sbid = find(sbid);
            prid = find(prid);
            
            if model.rev(Vex(irxn))
                lb = zeros(2,1); %kcat lb
                ub = zeros(2,1); %kcat ub     
                
                x0 =zeros(2,1); %kcat initial value
                
                lb_p = zeros(length(mc),1);
                ub_p = zeros(length(mc),1);
                
                x0_p = zeros(length(mc),1);
                x0_p(prid(mc(prid)>0)) = 1e-3;
                x0_p = x0_p(prid);
            else
                lb = zeros(1,1); %kcat lb
                ub = zeros(1,1); %kcat ub   
                
                x0 = zeros(1,1);
                
                lb_p = [];
                ub_p = [];
                
                x0_p = [];
            end
            ub(ub==0) = Inf;
            
            lb_s = zeros(length(mc),1);
            ub_s = zeros(length(mc),1);
            
            lb_s(sbid(mc(sbid)>0)) = 1e-9;
            lb_s(sbid(mc(sbid)==0)) = 1;
            lb_s = lb_s(sbid);
            ub_s(sbid(mc(sbid)>0)) = mc(sbid(mc(sbid)>0));  
            ub_s(sbid(mc(sbid)==0)) = 1;
            ub_s = ub_s(sbid);
            
            x0_s = zeros(length(mc),1);
            x0_s(sbid(mc(sbid)>0)) = 1e-3;
            x0_s = x0_s(sbid);
            
            if ~isempty(lb_p)
                lb_p(prid(mc(prid)>0)) = 1e-9;
                lb_p(prid(mc(prid)==0)) = 1;
                lb_p = lb_p(prid);
            end
            
            if ~isempty(ub_p)
                ub_p(prid(mc(prid)>0)) = mc(prid(mc(prid)>0));
                ub_p(prid(mc(prid)==0)) = 1;
                ub_p = ub_p(prid);
            end           
          
            lb = [lb;lb_s;lb_p];
            ub = [ub;ub_s;ub_p];           
            
            x0(x0==0) = 1;
            x0 = [x0;x0_s;x0_p];
            
%             lb = [0;0;1e-9;1e-9];
%             ub = [Inf;Inf;mc(sbid);mc(prid)];
%             x0 = [1;1;1e-3;1e-3];
            fprintf('%s\n',model.rxns{Vex(irxn)});            
            [x,fval,flag,output,lambda] =...
            fmincon(@vflux,x0,[],[],[],[],lb,ub,[]);
            if flag > 0
                if model.rev(Vex(irxn))
                    pvec.kcat_fwd(Vex(irxn)) = x(1);
                    pvec.kcat_bkw(Vex(irxn)) = x(2);
                    pvec.K(sbid,Vex(irxn)) = x(3);
                    pvec.K(prid,Vex(irxn)) = x(4);
                elseif ~model.rev(Vex(irxn))
                    pvec.kcat_fwd(Vex(irxn)) = x(1);
                    pvec.K(sbid,Vex(irxn)) = x(2);
                end
                pvec.Vmax(Vex(irxn)) = 1;
                flux = TKinetics(model,pvec,mc,Vex(irxn));
                fprintf('Vss = %3.6g \t flux = %3.6g\n',model.Vss(Vex(irxn)),flux);
            else
                fprintf('Infeasible Solution to NonLinear Optimization for %s\n',model.rxns{Vex(irxn)});
                fprintf('Continuing to next reaction\n');
                continue
            end
        end
    end
end

function ss_flux = vflux(x)
kcatfwd = x(1);
pvec_1 = pvec;

if model.rev(Vex(irxn))    
    kcatbkw = x(2);
    Ksub = x(3);
    Kprd = x(4);
    pvec_1.K(prid,Vex(irxn)) = Kprd;
    pvec_1.kcat_bkw(Vex(irxn)) = kcatbkw;
%     pr = mc(prid)/Kprd;
else
    Ksub = x(2);     
end

pvec_1.kcat_fwd(Vex(irxn)) = kcatfwd;

pvec_1.K(sbid,Vex(irxn)) = Ksub;

[~,vflux] = TKinetics(model,pvec_1,mc,Vex(irxn));

% 
% sr = mc(sbid)/Ksub;
% 
% if model.rev(Vex(irxn))
%     vflux = (kcatfwd*sr-kcatbkw*pr)./...
%         (1+sr+pr);
% elseif ~model.rev(Vex(irxn))
%     vflux = (kcatfwd*sr)./...
%         (1+sr);
% end

ss_flux = (vflux-model.Vss(Vex(irxn)))^2;


end
end

    