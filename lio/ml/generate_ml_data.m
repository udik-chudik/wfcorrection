N_TRAIN = 500;
N_VALID = 500;

XTdata = zeros(41,41,1,N_TRAIN);
YTdata = zeros(N_TRAIN, 1);
for i=1:N_TRAIN
    y = (0.5-rand(1, 1))*1000;
    YTdata(i) = y;
    x = takeImage([0 0], [0 0 0 y 0 0 0 0 0]*1e-9, 'IRS_180' , 0, 0);
    XTdata(:,:,1,i) = log10(x(128-20:128+20, 128-20:128+20));
    disp(i)
end


XVdata = zeros(41,41,1,N_VALID);
YVdata = zeros(N_VALID, 1);
for i=1:N_VALID
    y = (0.5-rand(1, 1))*1000;
    YVdata(i) = y;
    x = takeImage([0 0], [0 0 0 y 0 0 0 0 0]*1e-9, 'IRS_180' , 0, 0);
    XVdata(:,:,1,i) = log10(x(128-20:128+20, 128-20:128+20));
    disp(i)
end



%% Learn



miniBatchSize  = 50;
validationFrequency = floor(numel(YTdata)/miniBatchSize);
options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',30, ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',20, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{XVdata,YVdata}, ...
    'ValidationFrequency',100, ...
    'Plots','training-progress', ...
    'Verbose',false);

net = trainNetwork(XTdata,YVdata,layers_57,options);

imagesc(XVdata(:,:,1,1))

