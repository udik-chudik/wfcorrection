global errs;
global last_x;

last_x = 0;

wavelength = 5e-7;

errs = [];

global N_ACT;

N_ACT = 12;

x0 = zeros(N_ACT*N_ACT, 1);

%global til_tilt_coeff;

%til_tilt_coeff = [0 0];

%til_tilt_coeff = fminsearch(@(t) sum(sum(takeImageWithPlanetIRS(x0, t))), [0 0]);

%I = takeImage(x0, til_tilt_coeff);
%imagesc(I);

fmin(x0)

%options = optimset('PlotFcns',@optimplotfval);
%x = fminsearch(@fmin, x0, options);
%return

%'MutationFcn', {@mutationadaptfeasible, 2, 1}
options = optimoptions('ga','FunctionTolerance',1e-10,'PlotFcns', {@gaplotbestindiv,@gaplotbestf,@gaplotexpectation,@gaplotrange}, 'MutationFcn', {@mutationadaptfeasible, 100, 100});
x = ga(@fmin, N_ACT*N_ACT, [], [], [],[], ones(1,N_ACT*N_ACT)*-160, ones(1,N_ACT*N_ACT)*160,[],options);    % 3*RMS
%x = fmin_adam(@fmin, x0, 1);

%options = optimoptions('patternsearch','MeshTolerance',0.1,'StepTolerance',0.1, 'PlotFcns', {@psplotbestf, @psplotbestx});
%x = patternsearch(@fmin, x0,[], [], [],[], ones(1,N_ACT*N_ACT)*-40*3, ones(1,N_ACT*N_ACT)*40*3, [], options);

%options = optimoptions('surrogateopt','PlotFcn', @optimplotfval);
%x = surrogateopt(@fmin,ones(1,N_ACT*N_ACT)*-100*1e-9, ones(1,N_ACT*N_ACT)*100*1e-9);


%options = optimoptions('simulannealbnd','FunctionTolerance',1e-10,'PlotFcns', {@saplotbestx,@saplotbestf,@saplotx,@saplotf});
%x = simulannealbnd(@fmin, x0, ones(N_ACT*N_ACT, 1)*-100, ones(N_ACT*N_ACT, 1)*100,options);    % 3*RMS


%options = optimoptions('particleswarm','PlotFcn', @pswplotbestf, 'MaxStallIterations', 100);
%x = particleswarm(@fmin,N_ACT*N_ACT,ones(N_ACT*N_ACT, 1)*-40*3, ones(N_ACT*N_ACT, 1)*40*3,options);

img1 = takeImage([0 0], x0, 'IRS_180', 1, 1);
img2 = takeImage([0 0], x*1e-9, 'IRS_180', 1, 1);

imagesc([log10(cutZone(img1)) log10(cutZone(img2))]);
imagesc([log10(img1) log10(img2)]);

%imagesc([log10(cutZone(img1)) log10(cutZone(img2))]);

%rectangle('Position', [110 100 20 20])
%rectangle('Position', [110+256 100 20 20]);

%s1 = sum(cutZone( takeImage(reshape(x,[N_ACT N_ACT]), til_tilt_coeff)), 'all');
%s0 = sum(cutZone( takeImage(x0, til_tilt_coeff)), 'all');

%s0/s1



function s = fmin(x)
    global errs;
    global last_x;
    
    % with regularization
    a = 1;   % Желаемый контраст в зоне
    b = 10;     % среднее по изображению ~ 1e-9
    c = 1e-7;    % Средний квадрат отклонения ~ RMS^2 -> 2.5e-17 для 5 нм
    image = takeImage([0 0], x*10*1e-9, 'IRS_180', 1, 1);
    %s = a*mean(cutZone(image), 'all') + b*mean(image, 'all'); + c*dot(x,x)/length(x);
    %s = a*mean(log10(image(:,1:127)), 'all');
    s = sum(image, 'all');
    %s = sum(cutZone(image), 'all');
    errs = [errs s];
end