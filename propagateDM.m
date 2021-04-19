function [phase, amp] = propagateDM(act)
wavelength = 5e-7;
    gridsize = 512;
    diam = 0.1;
    focal_length = 0.3;
    beam_ratio = 0.2;

    wf = prop_begin( diam, wavelength, gridsize, beam_ratio );
    wf = prop_circular_aperture( wf, diam/2 );
    wf = prop_define_entrance( wf );

    tt = [-0.4962   -0.4981];

    s = size(act);
    tiptilt = zeros(512,512);
    for k=1:gridsize
        tiptilt(:,k) = linspace(0,wavelength*tt(1),gridsize);
    end
    
    for k=1:gridsize
        tiptilt(k,:) = tiptilt(k,:) + linspace(0,wavelength*tt(2),gridsize);
    end
    
    wf = prop_add_phase(wf, tiptilt);
    
    wf = prop_dm(wf, act, (s(1)-1)/2, (s(2)-1)/2, diam/s(1));

    
    wf = prop_propagate( wf, focal_length*100 );
    wf = prop_lens( wf, focal_length );
    wf = prop_propagate( wf, focal_length );


    
    
    
    
    amp = prop_get_amplitude(wf);
    phase = prop_get_phase(wf);
    
    amp1 = amp;
    amp2 = flip(amp);
    phase1 = phase;
    phase2 = flip(phase);
    
    E1 = amp1.*exp(1i*phase1);
    E2 = amp2.*exp(1i*(phase2+pi));
    
    E = E1 + E2;
    
    amp = abs(E);
    phase = angle(E);
    
    [wf, sampling] = prop_end( wf );
end

