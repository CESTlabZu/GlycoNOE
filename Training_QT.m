%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GLARES QUANTIFICATION NETWORK TRAINING 
%
% Authors: Leqi Yin, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% load training data here first

load('TMP_4p7T_0p25P_Je18-M.mat')% load sample testing data
%%
power=1;

k_cut=45:85;

Omax=1000;
step=25;
offset= -Omax:step:Omax;
k_0=[-2000, -1750, -1500, -1250, offset, 1250, 1500,1750,2000]';
gauss = 100;

randidx = 1:262144;%262144
num_samples = length(randidx);




%% Training QT

DataTrain1=cat(1,reshape(DL_CESTZ(k_cut,:),1,length(k_cut),[]),...
    reshape(DL_LoD_inv(k_cut,:),1,length(k_cut),[]));

DataOutput1=cat(1,reshape(DL_CESTZ_ref(k_cut,randidx),1,length(k_cut),[]),...
    reshape(DL_CESTZ_effMT(k_cut,randidx),1,length(k_cut),[]));

% DataTrain1=reshape(DL_CESTZ(k_cut,:),1,length(k_cut),[]);
% DataOutput1=reshape(DL_CESTZ_ref(k_cut,randidx),1,length(k_cut),[]);
% uncomment for trainning NN for ablation study

Train_mu=mean(DataTrain1(1,:,:),3);
Train_sigma=std(DataTrain1(1,:,:),[],3);

DataTrain=(DataTrain1-Train_mu)./Train_sigma;
DataOutput=(DataOutput1-Train_mu)./Train_sigma;

trainRatio = 0.8;
numTrain = floor(trainRatio * length(randidx));

XTrain1 = DataTrain(:,:,1:numTrain);
YTrain1 = DataOutput(:,:,1:numTrain);

XVal1 = DataTrain(:,:,numTrain+1:end);
YVal1 = DataOutput(:,:,numTrain+1:end);

XTraincell1 = squeeze(num2cell(XTrain1, [1 2]));
YTraincell1 = squeeze(num2cell(YTrain1, [1 2]));
XValcell1 = squeeze(num2cell(XVal1, [1 2]));
YValcell1 = squeeze(num2cell(YVal1, [1 2]));



XTraincell1 = cellfun(@transpose, XTraincell1, 'UniformOutput', false);
YTraincell1 = cellfun(@transpose, YTraincell1, 'UniformOutput', false);
XValcell1 = cellfun(@transpose, XValcell1, 'UniformOutput', false);
YValcell1 = cellfun(@transpose, YValcell1, 'UniformOutput', false);


%%

minLength = 41;
numFeatures = 2;      % One feature per time step
filterSize = 3;       % Kernel size
numFilters = 64;     % Number of filters in Conv layers

layers=GLARES_NN_QT(minLength, filterSize, numFilters);

%%

options = trainingOptions('adam', ...
    'MaxEpochs', 100, ...
    'MiniBatchSize', 16, ...
    'InitialLearnRate', 1e-4, ...
    'Shuffle', 'every-epoch', ...
    'Verbose', true, ...
    'ValidationData', {XValcell1, YValcell1}, ...
    'ValidationFrequency',5000,...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.3, ...
    'LearnRateDropPeriod', 25);

net_QT=trainnet(XTraincell1, YTraincell1, layers, "mse", options);

%%

% Example noisy input
sampleIdx = 615;
sampleIdx = randperm(size(XVal1,3),1);
sample = XValcell1{sampleIdx};


% Predict reference
result1 = predict(net_QT, sample);



% Plot comparison
figure;
plot(squeeze(XVal1(1,:,sampleIdx)).*Train_sigma+Train_mu, 'r--', 'LineWidth', 1); hold on;
plot(squeeze(YVal1(1,:,sampleIdx)).*Train_sigma+Train_mu, 'k-', 'LineWidth', 1);
plot(result1(:,1)'.*Train_sigma+Train_mu, 'b-', 'LineWidth', 1);


%% Testing

Data_Test1=squeeze(TMP_Zspectra_matrix(:,:,:,power));
Data_Test=B0correct_CEST(Data_Test1,32);
Data_Test_1D=reshape(Data_Test,32*32,89);

Test_mu=mean(Data_Test_1D,1);
Test_sigma=std(Data_Test_1D,[],1);

%%  LD


beta0_2pool= [0.9, 0, 280, 0.1, 0, 5000]; % initial test
lb_2pool=[0.02, -200, 60, 0, -800, 2000]; % lower bound
ub_2pool=[1, 200,2000,1, 800, 20000]; % upper bound

x=k_0;
x_2pool=[-4000
    -3500
    -3000
    -2500

    -200
    -150
    -100
    -50
    0
    50
    100
    150
    200

    2500
    3000
    3500
    4000]/2;
x_2pool=x_2pool';

for ii=1:32
    for jj=1:32
        sig1=double(squeeze(Data_Test(ii,jj,:)));
        sig=1-sig1';
        for i=1:1:4
            sig_2pool(i)=sig(i);
        end
        for i=5:1:13
            sig_2pool(i)=sig(i+36);
        end
        for i=14:1:17
            sig_2pool(i)=sig(i+72);
        end

        Delta=[1];

        options=optimset('lsqcurvefit') ;
        options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

        [beta_2pool,resnorm,residual,exitflag,output,lambda,jacobian] = ...
            lsqcurvefit(@matsolv_2pool, beta0_2pool, x_2pool, sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

        LoD_2pool_Test(ii,jj,:)=matsolv_LD(beta_2pool,x,Delta);

        sprintf("LoD----------------------- %d",ii)
    end
end

LoD_inv_Test=1-LoD_2pool_Test;
%%

Test_QT1=cat(3,reshape(Data_Test(:,:,k_cut),32,32,1,41),reshape(LoD_inv_Test(:,:,k_cut),32,32,1,41));

for i = 1:32
    for j =1:32
        Test_QT=(squeeze(Test_QT1(i,j,:,:))-Test_mu(k_cut))./Test_sigma(k_cut);
        TMP_pred_QT1=predict(net_QT, Test_QT');
        TMP_pred_QT(i,j,:)=TMP_pred_QT1(:,1)'.*Test_sigma(k_cut)+Test_mu(k_cut);
    end
end




%%

for ii=1:32
    for jj=1:32

        pred_arex(ii,jj,:) = ((1./Data_Test(ii,jj,k_cut)) - (1./TMP_pred_QT(ii,jj,:))).*TMP_R1W_cal_matrix(ii,jj)*(1+TMP_fm_matrix(ii,jj));

        pred_mtr(ii,jj,:) = (1-Data_Test(ii,jj,k_cut)) - (1-TMP_pred_QT(ii,jj,:,1));

    end
end

int_range=[50:56]-k_cut(1)+1;


for i=1:32
    for j=1:32
        xx = k_cut;
        pred_arex_area(i,j) = trapz(xx(int_range), pred_arex(i,j,int_range));
        pred_mtr_area(i,j) = trapz(xx(int_range), pred_mtr(i,j,int_range));
    end
end

%%
loss_N = abs(pred_arex_area-TMP_AREX_GT_int_matrix(:,:,power))./TMP_AREX_GT_int_matrix(:,:,power);
mean(loss_N(:))

loss_N_MTR = abs(pred_mtr_area-TMP_MTR_GT_int_matrix(:,:,power))./TMP_MTR_GT_int_matrix(:,:,power);
mean(loss_N_MTR(:))



%%

figure(1),imagesc(pred_arex_area,[0 0.3])

figure(2),imagesc(pred_mtr_area,[0 0.3])

figure(3),imagesc(TMP_AREX_GT_int_matrix(:,:,1),[0 0.3])



