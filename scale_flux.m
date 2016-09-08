function flux = scale_flux(flux)
for irxn = 1:size(flux,1)
    vf = flux(irxn,:);
    vf(abs(vf)<=1e-10) = 0;
    flux(irxn,:) = vf;
%     if vf>0 | vf<0
%         if abs(vf)<=1e-10
%             flux(irxn) = 0;
%         end    
%     end
end
            
            
        
           