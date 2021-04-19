C0 = zeros(3,3);

[I0, A0, P0] = takeImage(C0);

G_amp = zeros(512*512, 9);
G_phase = zeros(512*512, 9);

for k=1:9
    C = C0;
    C(imod(k,3), idiv(k,3)) = 5e-8;
    %disp([imod(k,3) idiv(k,3)])
    [img, amp, phase] = takeImage(C);
    amp = amp - A0;
    phase = phase - P0;
    
    for m=1:512*512
        G_amp(m, k) = amp(imod(m,512), idiv(m,512));
        G_phase(m, k) = phase(imod(m,512), idiv(m,512));
    end
    
end

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