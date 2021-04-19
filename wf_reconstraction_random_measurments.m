wavelength = 5e-7;

C0 = zeros(5,5);
I0 = takeImage(C0);

dls = [];
M = [];

for N=1:10
    Cp = generateRandomDM(C0)*wavelength/10;
    Ip = takeImage(Cp);
    Im = takeImage(-Cp);
    dls = [dls (Ip - Im)/2];
    
    p_abs = abs(sqrt((Ip+Im)/2 - I0));
    p_img = propagateDM(Cp);
    p_complex = p_abs.*exp(1i*p_img);
    M = [M; -imag(p_complex) real(p_complex) ];
    
end

re = zeros(512, 512);
im = zeros(512, 512);
for a=1:512
    for b=1:512
        dl = [];
        to_inv = [];
        for k=0:(N-1)
            dl = [dl dls(a, b + 512*k)];
            to_inv = [to_inv; M(a+512*k, b) M(a+512*k, b+512)];
        end
        c = pinv(to_inv)*dl';
        c = c/2;
        re(a,b) = c(1);
        im(a,b) = c(2);
    end
end


I_rec = (abs(re+1i*im).^0.2);

%{
for a=1:512
    for b=1:512
        if I_rec(a,b) > 0.05
            I_rec(a,b) = 0;
        end
    end
end
%}
imagesc(I_rec);
%imagesc(I_rec);
%imagesc(I0-abs(re+1i*im).^2);
%imagesc(angle(re+1i*im));