function flux = scale_flux(flux)
for irxn = 1:length(flux)
    vf = flux(irxn);
    if ~isnan(vf) && vf
        if abs(vf)<=1e-8
            flux(irxn) = 0;
        end    
    end
end
            
            
        
           