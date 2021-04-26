global errs;
wavelength = 5e-7;

errs = [];

global N_ACT;

N_ACT = 100;

x0 = zeros(N_ACT*N_ACT,1);

%x = fminsearch(@fmin, x0);
options = optimoptions('ga','FunctionTolerance',1e-10,'PlotFcn', @gaplotbestf);
x = ga(@fmin, N_ACT*N_ACT, [], [], [],[], ones(1,N_ACT*N_ACT)*-0.05*wavelength, ones(1,N_ACT*N_ACT)*0.05*wavelength,[],options);

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
    
    s = sum(sum(cutZone(takeImage(x))));
    
    errs = [errs s];
end