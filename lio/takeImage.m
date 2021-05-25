function [Ifinal, sampling] = takeImage(tt, x, coro_type, use_planet, use_errors)
% tt - [tip_tilt_x, tip_tilt_y]
% x = [....]
    N_ACT = sqrt(length(x));
    if ~(floor(N_ACT) == N_ACT) || N_ACT < 3
        error('---> takeImage: length(x) has to be square number and greater than 8');
    end
        


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
  
  wavefront   = telescope_with_dms(wavefront,f_lens, use_errors, reshape(x, [N_ACT N_ACT]));
  %wavefront   = telescope_with_dms_auto(wavefront,f_lens, use_errors, reshape(x, [N_ACT N_ACT]));
  
  switch coro_type
      case 'LIO'
          wavefront   = coronagraph_lio(wavefront, f_lens, 'GAUSSIAN', diam);
          [Ifinal, sampling] = prop_end(wavefront);
      case 'THRU'
          wavefront  = prop_lens(wavefront, f_lens );  % 'reimaging lens'
          wavefront  = prop_propagate(wavefront, f_lens ); % 'final focus'
          [Ifinal, sampling] = prop_end(wavefront);
      case 'IRS_180'
          Ifinal = coronagraph_rot(wavefront, f_lens);
      otherwise
          error('---> takeImage: unknown coronagraph type');
  end
  Ifinal = Ifinal / 0.071296289442758;  % Normalize intensity (THRU + without errors)
  
  if use_planet
     Iplanet = takeImage([-5*wavelength/diam 0], x, coro_type, 0, use_errors);
     Ifinal = Ifinal + Iplanet/1e4;
  end
  
end

