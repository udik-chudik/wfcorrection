function [wavefront, sampling] = takeImage(x)
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

  wavefront   = telescope_with_dms(wavefront,f_lens, 1, reshape(x, [49 49]));
  
  wavefront   = coronagraph(wavefront, f_lens, 'GAUSSIAN', diam);

  [wavefront, sampling] = prop_end(wavefront);
  
end

