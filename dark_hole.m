wavelength = 5e-7;
C1p = [wavelength/10 0 0; 0 wavelength/10 0; 0 0 wavelength/10];
C1m = -C1p;

C2p = [0 wavelength/10 wavelength/10; wavelength/10 0 wavelength/10; wavelength/10 wavelength/10 0];
C2m = -C2p;

C0 = zeros(3,3);

I0 = takeImage(C0);
I1p = takeImage(C1p);
I1m = takeImage(C1m);
I2p = takeImage(C2p);
I2m = takeImage(C2m);


d1 = (I1p - I1m)/2;
d2 = (I2p - I2m)/2;

p1_abs = abs(sqrt((I1p+I1m)/2 - I0));
p2_abs = abs(sqrt((I2p+I2m)/2 - I0));

p1_img = propagateDM(C1p);
p2_img = propagateDM(C2p);

dp1_complex = p1_abs.*exp(1i*p1_img);
dp2_complex = p2_abs.*exp(1i*p2_img);


%s = size(I0);
s = [41 41];

re = zeros(s);
im = zeros(s);
for a=1:s(1)
    for b=1:s(2)
        c = pinv([-imag(dp1_complex(a,b)) real(dp1_complex(a,b));-imag(dp2_complex(a,b)) real(dp2_complex(a,b))])*[d1(a,b) d2(a,b) ]';
        c = c/2;
        re(a,b) = c(1);
        im(a,b) = c(2);
    end
end

E_rec = re+1i*im;

%% Теперь собственно коррекция

N_act = 10;

C = zeros(N_act, N_act);
for jj=1:5

[Pm0, Am0] = propagateDM(C);

[img, A0, P0] = takeImage(C0);

P0 = cutZone(P0);
A0 = cutZone(A0);

acc0 = reshape(C,[],1);

G_complex = zeros(s(1)*s(2), N_act*N_act);

for k=1:N_act*N_act
    
    acc = acc0;
    acc(k) = acc(k) + wavelength/100;
    [phase, amp] = propagateDM(reshape(acc, [N_act,N_act]));
    
    phase = cutZone(phase);
    amp = cutZone(amp);
    
    %G_amp(:, k) = reshape(amp, [], 1);
    %G_phase(:, k) = reshape(phase, [], 1);
    
    Gc = amp.*exp(1i*phase) - cutZone(Am0).*exp(1i*cutZone(Pm0));
    G_complex(:, k) = reshape(Gc, [], 1);
end


%{
for k=1:9
    C = C0;
    C(imod(k,3), idiv(k,3)) = 1e-8;
    
    [phase, amp] = propagateDM(C);
    amp = cutZone(amp);
    phase = cutZone(phase);
    
    %amp = amp;
    %phase = phase;
    
    for m=1:s(1)*s(2)
        G_amp(m, k) = amp(imod(m,s(1)), idiv(m,s(2)));
        G_phase(m, k) = phase(imod(m,s(1)), idiv(m,s(2)));
    end
    
end
%}

%G0_complex = A0.*exp(1i*P0);
%G0_complex = reshape(G0_complex, [], 1);

%G_complex = (G_amp).*exp(1i*(G_phase));

G_inv = pinv(G_complex);
a_complex = G_inv*reshape(1i*E_rec, [], 1);


% Другая оценка
E_rec = A0.*exp(1i*P0);
a_complex = G_inv*reshape(1i*E_rec, [], 1);
% Конец другой оценки

C = reshape(real(a_complex), [N_act,N_act])/10;

disp(sum(sum(cutZone(takeImage(C*wavelength/100)))));
end

figure(1);
imagesc(cutZone(takeImage(C0)));
figure(2);
imagesc(cutZone(takeImage(C*wavelength/1000)));


%{
figure(1);
imagesc(takeImage(C0).^.2);
figure(2);
imagesc(takeImage(C*wavelength/1000).^.2);
%}

%imagesc((takeImage(C0)-takeImage(C*5e-8*20)));


%Ig = reshape(G_amp*reshape(C, [], 1), [512,512]);


function [a] = imod(c,d)
    a = mod(c,d);
    if (a == 0)
        a = d;
    end
end

function [a] = idiv(c,d)
    a = 1+(c - imod(c,d))/d;
end



















