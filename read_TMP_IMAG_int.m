clear; clc; close all;

load('TMP_4p7T_0p25P_JE18-M.mat')

%%
power=1;

Omax=2000;
step=50;

offset= -Omax:step:Omax;

k_offset=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000]*0.5;

if power == 1
    tt=0.25;
else
    tt=0.5;
end

Data_Test_O=squeeze(TMP_Zspectra_matrix(:,:,:,power));

%add noise
for ii =1:32
    for jj=1:32
        noiseFactor=1/200;
        Data_Test_N(ii, jj, :)=Data_Test_O(ii, jj, :)+noiseFactor * randn(1,1,89);
    end
end

%% B0 correction
Data_Test_CB0=B0correct_CEST(Data_Test_N,32);

%% Denoise
for i = 1:32
    for j = 1:32
        
        % Initialize selection
        selected_CESTZ = zeros(5,89);
        
        % Include the center element
        selected_CESTZ(1,:) = squeeze(Data_Test_CB0(i, j, :));
        
        % Compute distances to all other points in the matrix
        distances = zeros(32, 32);
        for m = 1:32
            for n = 1:32
                distances(m, n) = sqrt((m - i)^2 + (n - j)^2);
            end
        end
        
        % Flatten and sort distances while keeping track of indices
        [sorted_distances, sorted_indices] = sort(distances(:));
        
        % Select the 4 closest valid points (excluding itself)
        for k = 2:5 % Start from 2 to exclude itself (first index is itself)
            index = sorted_indices(k);
            [ni, nj] = ind2sub([32, 32], index);
            
            selected_CESTZ(k, :) = squeeze(Data_Test_CB0(ni, nj, :));

        end

        Data_Test_sets_matrix(i,j,:,:)=selected_CESTZ;
    end
end

DN_1D=reshape(Data_Test_CB0,32*32,89);
DN_mu=mean(DN_1D);
DN_sigma=std(DN_1D);

load('fitnet_2Ddenoise.mat')

for i = 1:32
    for j =1:32
        Test_denoise=(squeeze(Data_Test_sets_matrix(i,j,:,:))-DN_mu)./DN_sigma;
        pred_denoise1(i,j,:,:)=predict(net_2Ddenoise, Test_denoise').*DN_sigma'+DN_mu';
    end
end
Data_Test_denoise = squeeze(pred_denoise1(:,:,:,1));

%%
% Data_Test_denoise=Data_Test_CB0;
%only use when use clean data

%% multiple-pool Lorentzian fit.

for i=1:32
     for j=1:32
        
        sig=(1-double(squeeze(Data_Test_denoise(i,j,:))));
        
        R1W_AREX=TMP_R1W_cal_matrix(i,j);
        fm_AREX=TMP_fm_matrix(i,j);
        
        x =k_offset';

        % initial test
        beta0 = [0.9, 0, 280,        0.025, -700, 120,       0.01, -400, 200,        0.001, 200, 100,        0.02, 600, 500,         0.1, 0, 5000];
        % lower bound
        lb =    [0.02, -200, 60,     0, -800, 80,            0, -600, 100,           0, 100, 0,              0, 400, 200,            0, -800, 2000];
        % upper bound
        ub =    [1, 200, 2000,       0.2, -600, 600,         0.2, -200, 1000,        0.2, 300, 500,          1, 900, 1000,           1, 800, 20000];


        Delta=[1];
        options=optimset('lsqcurvefit') ;
        options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

        [beta,resnorm,residual,exitflag,output,lambda,jacobian] = ...
            lsqcurvefit(@matsolv, beta0, x, sig, lb, ub, options, Delta) ;


        % GCG
        beta_GCG=beta;
        sig_simur_GCG=matsolv(beta_GCG,x,Delta);
        beta_GCG(10)=0;
        sig_simur_ref_GCG=matsolv(beta_GCG,x,Delta);



        MPLF_AREX_GCG(i,j,:)=(1./(1-sig_simur_GCG)-1./(1-sig_simur_ref_GCG))*R1W_AREX*(1+fm_AREX);%use this one!!!


       
        sprintf("MPLF----------------------- %d",i)
            
     end
end

clearvars beta_GCG beta


%%

% multiple-pool Voigt fit.

for i=1:32
     for j=1:32
        
        sig=(1-double(squeeze(Data_Test_denoise(i,j,:))));
        
        R1W_AREX=TMP_R1W_cal_matrix(i,j);
        fm_AREX=TMP_fm_matrix(i,j);
        
        x =k_offset';

        % initial test
        beta0 = [0.9, 0, 280,        0.025, -700, 120,       0.01, -400, 200,        0.001, 200, 100,        0.02, 600, 500,         0.1, 0, 5000];
        % lower bound
        lb =    [0.02, -200, 60,     0, -800, 80,            0, -600, 100,           0, 100, 0,              0, 400, 200,            0, -800, 2000];
        % upper bound
        ub =    [1, 200, 2000,       0.2, -600, 600,         0.2, -200, 1000,        0.2, 300, 500,          1, 900, 1000,           1, 800, 20000];


        Delta=[1];
        options=optimset('lsqcurvefit') ;
        options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;
        
        [beta,resnorm,residual,exitflag,output,lambda,jacobian] = ...
            lsqcurvefit(@matsolv_Voigt, beta0, x, sig, lb, ub, options, Delta) ;


        % GCG
        beta_GCG=beta;
        sig_simur_GCG=matsolv_Voigt(beta_GCG,x,Delta);
        beta_GCG(10)=0;
        sig_simur_ref_GCG=matsolv_Voigt(beta_GCG,x,Delta);



        MPVF_AREX_GCG(i,j,:)=(1./(1-sig_simur_GCG)-1./(1-sig_simur_ref_GCG))*R1W_AREX*(1+fm_AREX);%use this one!!!


       
        sprintf("MPVF----------------------- %d",i)
            
     end
end

clearvars beta_GCG beta

%% Lorentzian Difference


beta0_2pool= [0.9, 0, 280, 0.1, 0, 5000]; % initial test
lb_2pool=[0.02, -200, 60, 0, -800, 2000]; % lower bound
ub_2pool=[1, 200,2000,1, 800, 20000]; % upper bound

x=k_offset;
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
        sig1=double(squeeze(Data_Test_denoise(ii,jj,:)));
        sig=1-sig1;
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
            lsqcurvefit(@matsolv_LD, beta0_2pool, x_2pool, sig_2pool, lb_2pool, ub_2pool, options, Delta) ;

        sig_simu_2pool=reshape(matsolv_LD(beta_2pool,x,Delta),1,1,89);

        LoD_2pool_bg(ii,jj,:)=1-sig_simu_2pool;

        LoD_2pool(ii,jj,:)=(1./Data_Test_denoise(ii,jj,:)-1./(1-sig_simu_2pool)).*TMP_R1W_cal_matrix(ii,jj)*(1+TMP_fm_matrix(ii,jj));
        LoD_2pool_mtr(ii,jj,:) =(1-Data_Test_denoise(ii,jj,:)) - sig_simu_2pool;

        sprintf("LD----------------------- %d",ii)
    end
end



%%

%Hybrid(LoD+MPLF)

% LoD_2pool(44:46,:)=0;
for i=1:32
     for j=1:32
                
        sig=double(squeeze(LoD_2pool(i,j,45:85)));
        
        
        x =k_offset(45:85)';

        % initial test
        beta0 = [0.001, 200, 100,        0.02, 600, 100];
        % lower bound
        lb =    [0, 100, 0,              0, 400, 50];
        % upper bound
        ub =    [0.2, 300, 500,          1, 900, 300];

        Delta=[1]; 
        options=optimset('lsqcurvefit') ; 
        options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;
        
        [beta,resnorm,residual,exitflag,output,lambda,jacobian] = ...
            lsqcurvefit(@matsolv_hybrid, beta0, x, sig, lb, ub, options, Delta) ;


        % GCG
        beta_GCG=beta;
        beta_matrix(i,j,:)=beta;
        sig_simur_GCG=matsolv_hybrid(beta_GCG,x,Delta);
        beta_GCG(1)=0;
        sig_simur_ref_GCG=matsolv_hybrid(beta_GCG,x,Delta);



        hybrid_fit(i,j,:)=sig_simur_GCG-sig_simur_ref_GCG;


       
        sprintf("HB----------------------- %d",i)
            
     end
end


%%

%PLOF and PVOF


k=k_offset./200;

FitParam.PeakOffset = 1; % -1 ppm for glycogen
FitParam.PeakRange = [k(50),k(56)];
FitParam.Magfield = 42.58*4.7;

FitParam.satpwr = tt; 
FitParam.tsat = 5; 


%PVOF
for i=1:32
    for j=1:32
    FitParam.R1 = double(TMP_R1W_cal_matrix(i,j));
    FitParam.fm = double(TMP_fm_matrix(i,j));
    tempdata=double(squeeze(Data_Test_denoise(i,j,:)));
    [fitresult,fitparam] = pvof(k,tempdata', FitParam);
    PVOF_fit(i,j,:) = fitresult.Arex;
    end
end


clearvars fitresult fitparam

%%

%PLOF
for i=1:32
    for j=1:32
    FitParam.R1 = double(TMP_R1W_cal_matrix(i,j));
    FitParam.fm = double(TMP_fm_matrix(i,j));
    tempdata=double(squeeze(Data_Test_denoise(i,j,:)));
    [fitresult,fitparam] = plof_3(k,tempdata', FitParam);
    PLOF_fit(i,j,:) = fitresult.Arex;
    end
end

%% DL

k_cut=45:85;
QT_1D=reshape(Data_Test_denoise,32*32,89);

QT_mu=mean(QT_1D,1);
QT_sigma=std(QT_1D,[],1);


load('fitnet_QT.mat')
Test_QT1=cat(3,reshape(Data_Test_denoise(:,:,k_cut),32,32,1,41),reshape(LoD_2pool_bg(:,:,k_cut),32,32,1,41));

for i = 1:32
    for j =1:32
        Test_QT=(squeeze(Test_QT1(i,j,:,:))-QT_mu(k_cut))./QT_sigma(k_cut);
        DL_pred_QT1=predict(net_QT, Test_QT');
        DL_pred_QT=reshape(DL_pred_QT1(:,1)'.*QT_sigma(k_cut)+QT_mu(k_cut),1,1,length(k_cut));
        DL_pred_QT_matrix(i,j,:)=DL_pred_QT;
        DL_pred_arex(i,j,:) = ((1./Data_Test_denoise(i,j,k_cut)) - (1./DL_pred_QT)).*TMP_R1W_cal_matrix(i,j)*(1+TMP_fm_matrix(i,j));
    end
end

%%

for i=1:32
    for j=1:32
DL_pred_arex_int(i,j) = trapz(1:7, DL_pred_arex(i,j,6:12));
LoD_2pool_int(i,j) = trapz(1:7, LoD_2pool(i,j,50:56));
MPLF_int(i,j) = trapz(1:7, MPLF_AREX_GCG(i,j,50:56));
MPVF_int(i,j) = trapz(1:7, MPVF_AREX_GCG(i,j,50:56));
hybrid_int(i,j) = trapz(1:7, hybrid_fit(i,j,6:12));
PLOF_int(i,j) = trapz(1:7, PLOF_fit(i,j,6:12));
PVOF_int(i,j) = trapz(1:7, PVOF_fit(i,j,6:12));

    end
end

   
GT_int=double(squeeze(TMP_AREX_GT_int_matrix(:,:,power)));




%%

DL_pred_arex_int_f = medfilt2(DL_pred_arex_int,[3 3]);
LoD_2pool_int_f    = medfilt2(LoD_2pool_int,[3 3]);
MPLF_int_f         = medfilt2(MPLF_int,[3 3]);
MPVF_int_f         = medfilt2(MPVF_int,[3 3]);
hybrid_int_f       = medfilt2(hybrid_int,[3 3]);
PLOF_int_f         = medfilt2(PLOF_int,[3 3]);
PVOF_int_f         = medfilt2(PVOF_int,[3 3]);


%%

figure
% cmin = min(GT_int(:));
% cmax = max(GT_int(:));
cmin=0;
cmax=0.18;

colormap('jet');
t = tiledlayout(1, 8, 'TileSpacing', 'none', 'Padding', 'compact');

ax1 = nexttile;
imagesc(DL_pred_arex_int_f);
axis image off;
caxis([cmin cmax]);

ax2 = nexttile;
imagesc(LoD_2pool_int_f);
axis image off;
caxis([cmin cmax]);

ax3 = nexttile;
imagesc(MPLF_int_f);
axis image off;
caxis([cmin cmax]);

ax4 = nexttile;
imagesc(MPVF_int_f);
axis image off;
caxis([cmin cmax]);

ax6 = nexttile;
imagesc(PLOF_int_f);
axis image off;
caxis([cmin cmax]);

ax7 = nexttile;
imagesc(PVOF_int_f);
axis image off;
caxis([cmin cmax]);

ax5 = nexttile;
imagesc(hybrid_int_f);
axis image off;
caxis([cmin cmax]);

ax8 = nexttile;
imagesc(GT_int);
axis image off;
caxis([cmin cmax]);


cb = colorbar;  % Adds a colorbar to the right of the tiled layout;
cb.Ticks = [0 0.06 0.12 0.18];
cb.TickLabels = {'0','0.06','0.12','0.18'};
cb.FontSize = 18;
set([ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8], 'FontSize', 18);
set(gcf, 'Color', 'w');


