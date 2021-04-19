
[x, fval] = fminsearch(@wopt, [0 0])

%{
val = zeros(21,21);

k = 1;
for a=0.4:0.01:0.6
    m = 1;
    for b = -0.1:0.01:0.1
        val(k,m) = max(max(takeImageTT(zeros(3,3), a, b)));
        m = m+1;
    end
    k = k+1;
end
%}

img = takeImageTT(zeros(3,3), x(1), x(2));
imagesc(img);

function answer = wopt(x)
    answer = max(max(takeImageTT(zeros(3,3), x(1), x(2))));
    %answer = sum(sum(takeImageTT(zeros(3,3), x(1), x(2))));
end

function [img, amp, phase] = takeImageTT(act, t1, t2)
    wavelength = 5e-7;
    gridsize = 512;
    diam = 0.1;
    focal_length = 0.3;
    beam_ratio = 0.2;

    wf = prop_begin( diam, wavelength, gridsize, beam_ratio );
    wf = prop_circular_aperture( wf, diam/2 );
    wf = prop_define_entrance( wf );

    tiptilt = zeros(512,512);
    for k=1:gridsize
        tiptilt(:,k) = linspace(0,wavelength*t1,gridsize);
    end
    
    for k=1:gridsize
        tiptilt(k,:) = tiptilt(k,:) + linspace(0,wavelength*t2,gridsize);
    end
    
    wf = prop_add_phase(wf, tiptilt);
    
    s = size(act);
    
    wf = prop_dm(wf, act, (s(1)-1)/2, (s(2)-1)/2, diam/s(1));
  
    rms_error = 1e-9;
    c_freq = 100.0;
    high_power = 11/3;
    wf = prop_psd_errormap(wf, rms_error, c_freq, high_power, 'rms');

    wf = prop_propagate( wf, focal_length*100 );
    
    wf = prop_lens( wf, focal_length );
    wf = prop_propagate( wf, focal_length );
    
    %imagesc((prop_get_amplitude(wf)).^0.5);
    
    img = prop_get_amplitude(wf).^2;
    amp = prop_get_amplitude(wf);
    phase = prop_get_phase(wf);
    
    amp1 = amp;
    amp2 = rot90(amp, 2);
    phase1 = phase;
    phase2 = rot90(phase, 2);
    
    E1 = amp1.*exp(1i*phase1);
    E2 = amp2.*exp(1i*(phase2+pi));
    
    E = E1 + E2;
    
    amp = abs(E);
    phase = angle(E);
    img = amp.^2;
    
    
    [wf, sampling] = prop_end( wf );
end

