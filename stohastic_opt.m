wavelength = 5e-7;
global errs;
global N_ACT;

errs = [];

N_ACT = 10;
x0 = zeros(N_ACT,N_ACT);
x = ga(@wopt, N_ACT*N_ACT, [], [], [],[], ones(1,N_ACT*N_ACT)*-0.1*wavelength, ones(1,N_ACT*N_ACT)*0.1*wavelength,[],1);
%[x,fval,exitflag,output] = fminsearch(@wopt, reshape(x0, [], 1));
%options = optimoptions('patternsearch','FunctionTolerance', 1e-9, 'MeshTolerance', 1e-9, 'StepTolerance', 1e-9);
%x = patternsearch(@wopt, reshape(x0, [], 1), [], [], [], [], ones(1,N_ACT*N_ACT)*-0.1*wavelength, ones(1,N_ACT*N_ACT)*0.1*wavelength, options);
%options = optimoptions('simulannealbnd','FunctionTolerance', 1e-9);
%x = simulannealbnd(@wopt, reshape(x0, [], 1),ones(1,N_ACT*N_ACT)*-1, ones(1,N_ACT*N_ACT)*1);


img1 = takeImage(x0);
img2 = takeImage(reshape(x, [N_ACT N_ACT]));

imagesc([log10(img1) log10(img2)]);

rectangle('Position', [110 100 20 20])
rectangle('Position', [110+256 100 20 20]);

s1 = sum(cutZone( takeImage(reshape(x,[N_ACT N_ACT]))), 'all');
s0 = sum(cutZone( takeImage(x0)), 'all');

s0/s1

function answer = wopt(x)
    global errs;
    global N_ACT;
    %x = x*0.1*5e-7;
    %answer = max(max(cutZone( takeImage(reshape(x,[10 10])))));
    answer = sum(cutZone( takeImage(reshape(x,[N_ACT N_ACT]))), 'all');
    errs = [errs answer];
    answer = answer + abs(sum(x))*1e5/length(x);
end