wavelength = 5e-7;
C1p = [wavelength/10 0 0; 0 wavelength/10 0; 0 0 wavelength/10]/20;
C1m = -C1p;

C2p = [0 wavelength/10 wavelength/10; wavelength/10 0 wavelength/10; wavelength/10 wavelength/10 0]/20;
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



%M = [-imag(dp1_complexl) real(dp1_complex);-imag(dp2_complexl) real(dp2_complex)];

re = zeros(512, 512);
im = zeros(512, 512);
for a=1:512
    for b=1:512
        %c = [-imag(dp1_complex(a,b)) real(dp1_complex(a,b));-imag(dp2_complex(a,b)) real(dp2_complex(a,b))]\[d1(a,b) d2(a,b)]';
        c = pinv([-imag(dp1_complex(a,b)) real(dp1_complex(a,b));-imag(dp2_complex(a,b)) real(dp2_complex(a,b))])*[d1(a,b) d2(a,b) ]';
        c = c/2;
        re(a,b) = c(1);
        im(a,b) = c(2);
    end
end

I_rec = (abs(re+1i*im).^2);

for a=1:512
    for b=1:512
        if I_rec(a,b) > 0.05
            I_rec(a,b) = 0;
        end
    end
end

imagesc(I_rec);
%imagesc(I_rec);
%imagesc(I0-abs(re+1i*im).^2);
%imagesc(angle(re+1i*im));