function wavefront = coronagraph_rot(wavefront, f_lens)
    

    wavefront  = prop_lens(wavefront, f_lens );  % 'coronagraph imaging lens'
    wavefront  = prop_propagate(wavefront, f_lens );  % 'snm', 'occulter'
    
    amp = prop_get_amplitude(wavefront)/2;
    ampr = rot90(amp);
    %ampr = flip(amp, 1);
    %phase = prop_get_phase(wavefront);
    %wf = prop_get_wavefront(wavefront);
    
    
    %wavefront = wf + rot90(wf+exp(1i*pi));
    wavefront = real(amp.^2 + ampr.^2 - 2*amp.*ampr);
    
    
    
end