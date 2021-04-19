function [phase] = backPropagate(act)
    wavelength = 5e-7;
    gridsize = 512;
    diam = 0.01;
    focal_length = 0.3;
    beam_ratio = 0.2;

    wf = prop_begin( diam, wavelength, gridsize, beam_ratio );
    
    wf = prop_circular_aperture( wf, diam/2 );
    wf = prop_define_entrance( wf );


    wf = prop_dm(wf, act, 1, 1, 0.002);

    wf = prop_lens( wf, focal_length );
    wf = prop_propagate( wf, focal_length );


    
    phase = prop_get_phase(wf);
    
    [wf, sampling] = prop_end( wf );
end

