function wavefront = coronagraph_rot_pupil(wavefront, f_lens)
    

    
    amp = prop_get_amplitude(wavefront)/2;
    ph = prop_get_phase(wavefront);
    
    ampr = circshift(rot90(amp, 2), [1, 1]);
    phr = circshift(rot90(ph, 2), [1, 1]);
    
    
    wavefront.wf = fftshift(    amp.*exp(1i*ph) + ampr.*exp(1i*(phr+pi)) );

%{

    amp = abs(wavefront.wf);
    ph = angle(wavefront.wf);
    
    ampr = circshift(rot90(amp, 2), [1, 1]);
    phr = circshift(rot90(ph, 2), [1, 1]);
    
    wavefront.wf = amp.*exp(1i*ph) + ampr.*exp(1i*(phr+pi));
%}
    
    wavefront  = prop_lens(wavefront, f_lens );  % 'coronagraph imaging lens'
    wavefront  = prop_propagate(wavefront, f_lens );  % 'snm', 'occulter'
    
    %wavefront  = prop_circular_aperture(wavefront, 10d0, 'norm');
    
    %amp = prop_get_amplitude(wavefront)/2;
    
    %ampr = rot90(amp);
    %ampr = circshift(rot90(amp, 2), [1, 1]);
    %ampr = flip(amp, 1);
    %phase = prop_get_phase(wavefront);
    %wf = prop_get_wavefront(wavefront);
    
    
    %wavefront = wf + rot90(wf+exp(1i*pi));
    %wavefront = real(amp.^2 + ampr.^2 - 2*amp.*ampr) + 1e-17;   % 1e-18 presision of the comp.
    
    wavefront = prop_get_amplitude(wavefront).^2 + 1e-17;
    
    
    
end