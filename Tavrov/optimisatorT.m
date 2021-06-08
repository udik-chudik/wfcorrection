global errs;
global last_x;

last_x = 0;

wavelength = 5e-7;

errs = [];

global N_ACT;

N_ACT = 34;

x0 = zeros(N_ACT*N_ACT, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tavrov begin
                                    % Four-step phase shifting technique
    xT=zeros(N_ACT, N_ACT);         % \/����� �� ���� ������ ������� ����
    phase_offset=pi/7;              % phase off-set [rad];
    xT_(1)=phase_offset+0;          % four-step ... 0
    xT_(2)=phase_offset+pi/2;       % four-step ... 90
    xT_(3)=phase_offset+pi;         % four-step ... 180
    xT_(4)=phase_offset+3/2*pi;     % four-step ... 270
        xT_RadToM= wavelength/(2*pi);              % coeff [rad] to [m]
%        x_HalfPup=zeros(N_ACT, floor(N_ACT/2));   % at half of pupil
         x_HalfPup=zeros(floor(N_ACT/2),N_ACT);   % at half of pupil
    for Ti=1:length(xT_);           % step Ti/4
        x_HalfPup_Ti=x_HalfPup+xT_(Ti);
        xT_Ti=xT;
       % xT_Ti(1:N_ACT, 1:floor(N_ACT/2))=x_HalfPup_Ti;
       xT_Ti( 1:floor(N_ACT/2), 1:N_ACT)=x_HalfPup_Ti;
        xT_T(Ti,:,:)=xT_Ti.*xT_RadToM; % ?? normalization to [m]
        disp(Ti)
    end                             % Ti

%TTT=reshape(xT_T(1,:,:),N_ACT*N_ACT, 1);
imgT0 = takeImage([0 0], x0, 'IRS_180', 1, 0);
imgT1 = takeImage([0 0], reshape(xT_T(1,:,:),N_ACT*N_ACT, 1), 'IRS_180', 1, 1);
imgT2 = takeImage([0 0], reshape(xT_T(2,:,:),N_ACT*N_ACT, 1), 'IRS_180', 1, 1);
imgT3 = takeImage([0 0], reshape(xT_T(3,:,:),N_ACT*N_ACT, 1), 'IRS_180', 1, 1);
imgT4 = takeImage([0 0], reshape(xT_T(4,:,:),N_ACT*N_ACT, 1), 'IRS_180', 1, 1);
figure(10), imagesc(log10([imgT0, imgT1, imgT2 imgT3, imgT4])), colorbar;

%>>>>function [Ifinal, sampling] = takeImage(tt, x, coro_type, use_planet, use_errors)
imgT0 = takeImage([0 0], x0, 'THRU', 1, 1);
imgT0C = takeImage([0 0], x0, 'IRS_180', 1, 1);
figure(20), imagesc(log10([imgT0, imgT0C])), colorbar;


global use_errors_custom; % custom WF errors
use_errors_custom.rms_error =  40.0d-11;    % 40.0d-9; RMS wavefront error
use_errors_custom.c_freq =  15.0d0;         % correlation frequency (cycles / m)
use_errors_custom.high_power = 3.0d0;       % high frequency falloff
use_errors_custom.flnm = 'telescope_40Tnm.fits'; % filename
imgT0 = takeImage([0 0], x0, 'THRU', 1, 2);
imgT0C = takeImage([0 0], x0, 'IRS_180', 1, 2);

figure(20), imagesc(log10([imgT0, imgT0C])), colorbar;


imgT1 = takeImage([0 0], reshape(xT_T(1,:,:),N_ACT*N_ACT, 1), 'THRU', 1, 1);
imgT2 = takeImage([0 0], reshape(xT_T(2,:,:),N_ACT*N_ACT, 1), 'THRU', 1, 1);
imgT3 = takeImage([0 0], reshape(xT_T(3,:,:),N_ACT*N_ACT, 1), 'THRU', 1, 1);
imgT4 = takeImage([0 0], reshape(xT_T(4,:,:),N_ACT*N_ACT, 1), 'THRU', 1, 1);

figure(10), imagesc(log10([imgT0, imgT1, imgT2 imgT3, imgT4])), colorbar;
figure(11), imagesc(reshape(xT_T(2,:,:),35,35)), colorbar
size(xT_T(2,:,:))
% image algebra
Tphase=atan((imgT4.^2-imgT2.^2)./(imgT1.^2-imgT3.^2));
 
Tphase(isnan(Tphase))=0; %(1-isnan(Tphase)).*Tphase;
figure(11), imagesc(Tphase), colorbar, title('Phase in image plane');


figure(12), imagesc(angle(fftshift(ifft2(imgT1.*exp(Tphase))))), colorbar,  title('Phase in pupil plane');


Tvisib=(2*sqrt((imgT4.^2-imgT2.^2).^2+(imgT1.^2-imgT3.^2).^2)./(imgT1.^2+imgT2.^2+...
    imgT3.^2+imgT4.^2);
%figure(12), 


%figure(10), imagesc([log10(imgT0) log10(imgT1)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tavrov end

%global til_tilt_coeff;

%til_tilt_coeff = [0 0];

%til_tilt_coeff = fminsearch(@(t) sum(sum(takeImageWithPlanetIRS(x0, t))), [0 0]);

%I = takeImage(x0, til_tilt_coeff);
%imagesc(I);



%options = optimset('PlotFcns',@optimplotfval);
%x = fminsearch(@fmin, x0, options);
%return

%'MutationFcn', {@mutationadaptfeasible, 2, 1}
options = optimoptions('ga','FunctionTolerance',1e-10,'PlotFcns', {@gaplotbestindiv,@gaplotbestf,@gaplotexpectation,@gaplotrange}, 'MutationFcn', {@mutationadaptfeasible, 100, 100});
x = ga(@fmin, N_ACT*N_ACT, [], [], [],[], ones(1,N_ACT*N_ACT)*-40*3, ones(1,N_ACT*N_ACT)*40*3,[],options);    % 3*RMS
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
    %s = sum(image, 'all');
    s = sum(cutZone(image), 'all');
    errs = [errs s];
end