clear;
clc;
tic

load("trainingdata_4p7T_0p25P_Ab_AW_JE4.mat")

%%

% signal-to-noise ratio
i_SNR = 5;

% offset cut
k_cut = 1:89; 

n_features = length(k_cut);
n_outputs = 2;

randidx = 1:262144;

% Target matrix
DataTrain = DL_CESTZ;
DataOutput = [DL_amp.*10 DL_width./10000];

trainRatio = 0.8;
numTrain = floor(trainRatio * length(randidx));

XTrain = DataTrain(1:numTrain,k_cut);
YTrain = DataOutput(1:numTrain,:);
XVal = DataTrain(numTrain+1:end,k_cut);
YVal = DataOutput(numTrain+1:end,:);

XTrainCells = arrayfun(@(i) XTrain(i,:)', (1:size(XTrain,1))', 'UniformOutput', false);
YTrainCells = arrayfun(@(i) YTrain(i,:)', (1:size(YTrain,1))', 'UniformOutput', false);
XValCells = arrayfun(@(i) XVal(i,:)', (1:size(XVal,1))', 'UniformOutput', false);
YValCells = arrayfun(@(i) YVal(i,:)', (1:size(YVal,1))', 'UniformOutput', false);

%%

fc_layers = [
    featureInputLayer(n_features,Normalization="zerocenter")

    fullyConnectedLayer(512,"WeightsInitializer","he")
    batchNormalizationLayer
    eluLayer
    dropoutLayer(0.25)

    fullyConnectedLayer(256,"WeightsInitializer","he")
    batchNormalizationLayer
    eluLayer
    dropoutLayer(0.15)

    fullyConnectedLayer(64,"WeightsInitializer","he")
    batchNormalizationLayer
    eluLayer

    fullyConnectedLayer(16,"WeightsInitializer","he")
    batchNormalizationLayer
    eluLayer

    fullyConnectedLayer(n_outputs) % Predict amplitude and width
];

% conv_layers = [...
%     sequenceInputLayer(n_features,Normalization="zerocenter")
% 
%     convolution1dLayer(3,64,"Padding","same","WeightsInitializer","he")
%     eluLayer
% 
%     convolution1dLayer(3,32,"Padding","same","WeightsInitializer","he")
%     eluLayer
% 
%     convolution1dLayer(3,16,"Padding","same","WeightsInitializer","he")
%     eluLayer
% 
%     flattenLayer
% 
%     dropoutLayer(0.15)
% 
%     fullyConnectedLayer(256)
% 
%     fullyConnectedLayer(128)
% 
%     fullyConnectedLayer(n_outputs)
% 
%     regressionLayer
% 
% 
% ];





options = trainingOptions('adam', ...
    'MaxEpochs', 100, ...  % Increased for better learning
    'MiniBatchSize', 16, ... % Adjusted batch size for better updates
    'InitialLearnRate', 1e-4, ...
    'Shuffle', 'every-epoch', ...
    'Verbose', true, ...
    'ValidationData', {XValCells, YValCells}, ...
    'ValidationFrequency',5000,...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.3, ...
    'LearnRateDropPeriod', 25);


% net_QT=trainnet(XTrain, YTrain, conv_layers, "mse", options);	
net_QT=trainNetwork(XTrainCells, YTrainCells, conv_layers, options);



% uncomment the following lines to enable predictions on the testing data

% if trainingdata == 2
%     % revert the prediction results to their original scale
%     DL_pred = predict(net, DL_Input_testing).*300;
% else
%     DL_pred = predict(net, DL_Input_testing);
% end



