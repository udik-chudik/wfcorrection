global errs;
global last_x;

last_x = 0;

wavelength = 5e-7;

errs = [];

global N_ACT;

N_ACT = 32;

x0 = zeros(N_ACT*N_ACT,1);

global til_tilt_coeff;

til_tilt_coeff = fminsearch(@(t) sum(sum(takeImage(x0, t))), [0 0]);


%x = fminsearch(@fmin, x0);
options = optimoptions('ga','FunctionTolerance',1e-10,'PlotFcn', @gaplotbestf);
x = ga(@fmin, N_ACT*N_ACT, [], [], [],[], ones(1,N_ACT*N_ACT)*-3*1000e-9, ones(1,N_ACT*N_ACT)*3*1000e-9,[],options);    % 3*RMS

img1 = takeImage(x0, til_tilt_coeff);
img2 = takeImage(reshape(x, [N_ACT N_ACT]), til_tilt_coeff);

imagesc([log10(img1) log10(img2)]);

rectangle('Position', [110 100 20 20])
rectangle('Position', [110+256 100 20 20]);

s1 = sum(cutZone( takeImage(reshape(x,[N_ACT N_ACT]), til_tilt_coeff)), 'all');
s0 = sum(cutZone( takeImage(x0, til_tilt_coeff)), 'all');

s0/s1


function s = fmin(x)
    global errs;
    global last_x;
    last_x = x;
    global til_tilt_coeff;
    %s = sum(sum(cutZone(takeImage(x))));
    % with regularization
    a = 1;   % Желаемый контраст в зоне
    b = 1;     % среднее по изображению ~ 1e-9
    c = 1e8;    % Средний квадрат отклонения ~ RMS^2 -> 2.5e-17 для 5 нм
    image = takeImage(x, til_tilt_coeff);
    s = a*mean(cutZone(image), 'all') + b*mean(image, 'all') + c*dot(x,x)/length(x);
    
    errs = [errs s];
end