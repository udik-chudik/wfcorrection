N_TRAIN = 5000;
N_VALID = 5000;

XTdata = zeros(41,41,1,N_TRAIN);
YTdata = zeros(N_TRAIN, 8);
for i=1:N_TRAIN
    y = rand(8, 1);
    YTdata(i,:) = y;
    x = takeImage([0 0], [y(1)*500 y(2)*500 y(3)*500 y(4)*500 0 y(5)*500 y(6)*500 y(7)*500 y(8)*500]*1e-9, 'IRS_180' , 0, 0);
    XTdata(:,:,1,i) = x(128-20:128+20, 128-20:128+20)./0.12;
    disp(i)
end


XVdata = zeros(41,41,1,N_VALID);
YVdata = zeros(N_VALID, 8);
for i=1:N_VALID
    y = rand(8, 1);
    YVdata(i,:) = y;
    x = takeImage([0 0], [y(1)*500 y(2)*500 y(3)*500 y(4)*500 0 y(5)*500 y(6)*500 y(7)*500 y(8)*500]*1e-9, 'IRS_180' , 0, 0);
    XVdata(:,:,1,i) = x(128-20:128+20, 128-20:128+20)./0.12;
    disp(i)
end



%% Learn



miniBatchSize  = 50;
validationFrequency = floor(numel(YTdata)/miniBatchSize);
%{
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
%}

options = trainingOptions('sgdm', ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.2, ...
    'LearnRateDropPeriod',5, ...
    'MaxEpochs',20, ...
    'MiniBatchSize',64, ...
    'InitialLearnRate', 0.0001, ...
    'ValidationData',{XVdata,YVdata},...
    'ValidationFrequency',50,...
    'Plots','training-progress');

net = trainNetwork(XTdata,YTdata,layers_41,options);

pd = predict(net, XVdata);
scatter(pd(:,4), YVdata(:,4));
