function flux = scale_flux(flux)
for irxn = 1:length(flux)
    vf = flux(irxn);
    if vf>0 | vf<0
        if abs(vf)<=1e-10
            flux(irxn) = 0;
        end    
    end
end
            
            
        
           