re = zeros(512, 512);
im = zeros(512, 512);
for a=1:512
    for b=1:512
        %[re(a,b) im(a,b)] = inv([-imag(p1(a,b)) real(p1(a,b));-imag(p2(a,b)) real(p2(a,b))])*[d1 d2];
        c = [-imag(p1(a,b)) real(p1(a,b));-imag(p2(a,b)) real(p2(a,b))]\[d1(a,b) d2(a,b)]';
        re(a,b) = c(1);
        im(a,b) = c(2);
    end
end