function [wavefront, sampling] = takeImage(x, tt)

  global N_ACT;
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  diam        =    0.1d0;               % diameter (m)
  f_lens      =   24.0d0 * diam;        % focal length (m)
  beam_ratio  =    0.3d0;               % beam diameter fraction
  wavelength = 5e-7;
  grid_size = 256;
  
  wavefront   = prop_begin(diam, wavelength, grid_size, beam_ratio);
  wavefront   = prop_circular_aperture(wavefront, diam / 2.0d0);
  wavefront   = prop_define_entrance(wavefront);

  wavefront = add_tip_tilt(wavefront, tt(1), tt(2));
  wavefront   = telescope_with_dms(wavefront,f_lens, 1, reshape(x, [N_ACT N_ACT]));
  
  %wavefront   = coronagraph_lio(wavefront, f_lens, 'GAUSSIAN', diam);
  
  wavefront   = coronagraph_rot(wavefront, f_lens);

  %[wavefront, sampling] = prop_end(wavefront);
  
  %wavefront = wavefront/4;
  %wavefront = abs(wavefront).^2;
  %wavefront = wavefront/4;
 
end

