clear;clc

%load clean data
load('TMP_4p7T_0p25P_Je18-M.mat')

%%


Omax=2000;
step=50;

tt=0.25;

offset= -Omax:step:Omax;

k_offset=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000]*0.5;

Data_Test_O=TMP_Zspectra_matrix;

%% B0 correction
Data_Test_CB0=B0correct_CEST(Data_Test_O,32);

%% PLOF
k=k_offset./200;

FitParam.PeakOffset = 1; % -1 ppm for glycogen
FitParam.PeakRange = [k(50),k(56)];
FitParam.Magfield = 42.58*4.7;

FitParam.satpwr = tt; 
% FitParam.tsat = 5; 

for i=1:32
    for j=1:32
    FitParam.R1 = double(TMP_R1W_cal_matrix(i,j));
    FitParam.fm = double(TMP_fm_matrix(i,j));
    tempdata=double(squeeze(Data_Test_CB0(i,j,:)));
    [fitresult,fitparam] = plof_3(k,tempdata', FitParam);
    coefficient1_matrix(i,j,:)=fitresult.Coefficents1;
    coefficient2_matrix(i,j,:)=fitresult.Coefficents2;
    PLOF_fit(i,j,:) = fitresult.Arex;
    end
end