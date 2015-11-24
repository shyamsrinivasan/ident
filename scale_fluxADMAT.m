function flux = scale_fluxADMAT(flux)
for irxn = 1:length(flux)
    vf = flux(irxn);
    if ~isnan(vf) && getval(vf)
        if abs(vf)<=1e-10
            flux(irxn) = 0;
        end    
    end
end
            
            
        
           