global errs;
global last_x;

last_x = 0;

wavelength = 5e-7;

errs = [];

global N_ACT;

N_ACT = 42;

x0 = zeros(N_ACT*N_ACT,1);

%x = fminsearch(@fmin, x0);
options = optimoptions('ga','FunctionTolerance',1e-10,'PlotFcn', @gaplotbestf);
x = ga(@fmin, N_ACT*N_ACT, [], [], [],[], ones(1,N_ACT*N_ACT)*-3*500e-9, ones(1,N_ACT*N_ACT)*3*500e-9,[],options);    % 3*RMS

img1 = takeImage(x0);
img2 = takeImage(reshape(x, [N_ACT N_ACT]));

imagesc([log10(img1) log10(img2)]);

rectangle('Position', [110 100 20 20])
rectangle('Position', [110+256 100 20 20]);

s1 = sum(cutZone( takeImage(reshape(x,[N_ACT N_ACT]))), 'all');
s0 = sum(cutZone( takeImage(x0)), 'all');

s0/s1


function s = fmin(x)
    global errs;
    global last_x;
    last_x = x;
    %s = sum(sum(cutZone(takeImage(x))));
    % with regularization
    a = 2000;
    b = 20;
    c = 1e6;
    
    image = takeImage(x);
    s = a*mean(cutZone(image), 'all') + b*mean(image, 'all') + c*dot(x,x);
    
    errs = [errs s];
end