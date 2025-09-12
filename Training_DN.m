%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GLARES DENOISING NETWORK TRAINING 
%
% Authors: Leqi Yin, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% load training data and testing data here first


%%

Omax=2000;
step=50;

offset= -Omax:step:Omax;
k_0=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000]*0.5;
%262144
randidx = 1:262144;
num_samples = length(randidx);



Train_mu=reshape(mean(squeeze(DL_CESTZ_noised_sets(1,:,:)),2),1,89,1);
Train_sigma=reshape(std(squeeze(DL_CESTZ_noised_sets(1,:,:)),[],2),1,89,1);


%%

DataTrain = (DL_CESTZ_noised_sets(:,:,randidx)-Train_mu)./Train_sigma;
DataOutput = (DL_CESTZ_sets(:,:,randidx)-Train_mu)./Train_sigma;

trainRatio = 0.8;
numTrain = floor(trainRatio * length(randidx));

XTrain = DataTrain(:,:,1:numTrain);
YTrain = DataOutput(:,:,1:numTrain);

XVal = DataTrain(:,:,numTrain+1:end);
YVal = DataOutput(:,:,numTrain+1:end);

XTraincell = squeeze(num2cell(XTrain, [1 2]));
YTraincell = squeeze(num2cell(YTrain, [1 2]));
XValcell = squeeze(num2cell(XVal, [1 2]));
YValcell = squeeze(num2cell(YVal, [1 2]));

XTraincell = cellfun(@transpose, XTraincell, 'UniformOutput', false);
YTraincell = cellfun(@transpose, YTraincell, 'UniformOutput', false);
XValcell = cellfun(@transpose, XValcell, 'UniformOutput', false);
YValcell = cellfun(@transpose, YValcell, 'UniformOutput', false);


%% 
minLength=89;
filterSize = 3;       % Kernel size
numFilters = 256;     % Number of filters in Conv layers


layers=GLARES_NN_DN(minLength, filterSize, numFilters);

%%

options = trainingOptions('adam', ...
    'MaxEpochs', 100, ...  % Increased for better learning
    'MiniBatchSize', 32, ... % Adjusted batch size for better updates
    'InitialLearnRate', 1e-4, ...
    'Shuffle', 'every-epoch', ...
    'Verbose', true, ...
    'ValidationData', {XValcell, YValcell}, ...
    'ValidationFrequency',5000,...
        'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.3, ...
    'LearnRateDropPeriod', 25);
    net_2Ddenoise=trainnet(XTraincell, YTraincell, layers, "mse", options);


%% Sample Test
sampleIdx = randperm(size(XVal,3),1);
sample = XValcell{sampleIdx};
result = predict(net_2Ddenoise, sample);
figure;
plot(squeeze(XVal(1,:,sampleIdx).*Train_sigma+Train_mu), 'r--', 'LineWidth', 1); hold on;
plot(squeeze(YVal(1,:,sampleIdx).*Train_sigma+Train_mu), 'b-', 'LineWidth', 1);
plot(squeeze(result(:,1)'.*Train_sigma+Train_mu), 'k-', 'LineWidth', 1);


%% 

count=0;

invivo_1D=reshape(TMP_invivo_CESTZ_CB0,128*128,89);
for nn = 1:length(invivo_1D)
    if any(invivo_1D(nn, 50) ~= 0)
        count=count+1;
        non_zero_invivo_1D(count,:) = invivo_1D(nn,:);
    end
end
invivo_mu=mean(non_zero_invivo_1D);
invivo_sigma=std(non_zero_invivo_1D,[],1);





%% TMP denoise test


for i = 1:128
    for j =1:128
        if any(squeeze(TMP_invivo_CESTZ(i,j,:)) ~= 0)
            Test_denoise=(squeeze(TMP_invivo_CESTZ_sets_CB0(i,j,:,:))-invivo_mu)./invivo_sigma;
            pred_denoise1(i,j,:,:)=predict(net_2Ddenoise, Test_denoise').*invivo_sigma'+invivo_mu';
        else
            pred_denoise1(i,j,:,:)=zeros(89,5);
        end

    end
end
TMP_pred_denoise = squeeze(pred_denoise1(:,:,:,1));
%% sample Test Denoise

figure;

rows = randi(128);
cols = randi(128);

TMP_CESTZ_GT=squeeze(TMP_invivo_CESTZ_ns(rows,cols,:));

plot(squeeze(TMP_pred_denoise(rows,cols,:)), 'r-', 'LineWidth', 1); hold on;
plot(TMP_CESTZ_GT, 'k-', 'LineWidth', 1);

%%
sum(mse(TMP_pred_denoise,TMP_invivo_CESTZ))
