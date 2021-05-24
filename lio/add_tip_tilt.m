function wavefront = add_tip_tilt(wavefront, a, b)

    npa = prop_get_gridsize();
    
    phase = zeros(npa, npa);
    for i=1:npa
        for j=1:npa
            phase(i,j) = (a/npa)*(j-1) + (b/npa)*(i-1);
        end
    end

    wavefront = prop_add_phase(wavefront, phase);
end