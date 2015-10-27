function pvec = getRKparameter(model,pvec,mc,Vred)

h2o = find(strcmpi(model.mets,'h2o[c]'));
pic = find(strcmpi(model.mets,'pi[c]'));
pie = find(strcmpi(model.mets,'pi[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
he = find(strcmpi(model.mets,'h[e]'));

for irxn = 1:length(Vred)
    if model.Vss(Vred(irxn)) ~=0
        sbid = logical(model.S(:,Vred(irxn))<0);
        prid = logical(model.S(:,Vred(irxn))>0);
        if any(sbid) && any(prid)
            sbid([pic pie hc he h2o]) = 0;
            prid([pic pie hc he h2o]) = 0;
            
            sbid = find(sbid);
            prid = find(prid);
            
            if model.rev(Vred(irxn))
                lb = zeros(3,1); %[kcat Vmax] lb
                ub = zeros(3,1); %[kcat Vmax] ub     
                
                x0 =zeros(3,1); %kcat initial value
                
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
            
            fprintf('%s\n',model.rxns{Vred(irxn)});            
            [x,fval,flag,output,lambda] =...
            fmincon(@vflux,x0,[],[],[],[],lb,ub,[]);
            display(x);
        end        
    end
end

function ss_flux = vflux(x)
pvec_1 = pvec;
pvec_1.kcat_fwd(Vred(irxn)) = x(1);

if model.rev(Vred(irxn))
    pvec_1.kcat_bkw(Vred(irxn)) = x(2);
    pvec_1.K(sbid,Vred(irxn)) = x(3);
    pvec_1.K(prid,Vred(irxn)) = x(4);
else
    pvec_1.K(sbid,Vred(irxn)) = x(2);
end

[~,vflux] = TKinetics(model,pvec_1,mc,Vred(irxn));

ss_flux = (vflux-model.Vss(Vred(irxn)))^2;
end    
end
