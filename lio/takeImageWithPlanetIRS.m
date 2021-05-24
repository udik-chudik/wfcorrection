function [Ifinal, sampling] = takeImageWithPlanetIRS(x, tt)

  global N_ACT;

  diam        =    0.1d0;               % diameter (m)
  f_lens      =   24.0d0 * diam;        % focal length (m)
  beam_ratio  =    0.3d0;               % beam diameter fraction
  wavelength = 5e-7;
  grid_size = 256;
  
  wavefront   = prop_begin(diam, wavelength, grid_size, beam_ratio);
  %planet_wf = ones(grid_size,grid_size).*10.*exp(1i*repmat(linspace(0,pi*2*80, grid_size), grid_size,1));
  
  %wavefront = prop_add_wavefront(wavefront, planet_wf);
  
  wavefront   = prop_circular_aperture(wavefront, diam / 2.0d0);
  wavefront   = prop_define_entrance(wavefront);
    wavefront = add_tip_tilt(wavefront, tt(1), tt(2));
  wavefront   = telescope_with_dms(wavefront,f_lens, 1, reshape(x, [N_ACT N_ACT]));
  
  Icoro   = coronagraph_rot(wavefront, f_lens);

  
  
  %[Icoro, sampling] = prop_end(wavefront);
  
  
  wavefront   = prop_begin(diam, wavelength, grid_size, beam_ratio);
  wavefront = prop_add_phase(wavefront,  linspace(0, 5*wavelength/diam, grid_size));
  %planet_wf = ones(grid_size,grid_size).*10.*exp(1i*repmat(linspace(0,pi*2*80, grid_size), grid_size,1));
  
  %wavefront = prop_add_wavefront(wavefront, planet_wf);
  
  wavefront   = prop_circular_aperture(wavefront, diam / 2.0d0);
  wavefront   = prop_define_entrance(wavefront);
 wavefront = add_tip_tilt(wavefront, tt(1), tt(2));
  wavefront   = telescope_with_dms(wavefront,f_lens, 1, reshape(x, [N_ACT N_ACT]));
  
  Iplanet   = coronagraph_rot(wavefront, f_lens);
 
  %[Iplanet, sampling] = prop_end(wavefront);
  
  In = 0.0158;%max(Icoro, [], 'all');
  
  Ifinal = Icoro + Iplanet*1e0/1e4;
  Ifinal = Ifinal/In;
  
end