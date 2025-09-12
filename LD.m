clear;clc

%load clean data
load('TMP_4p7T_0p25P_Je18-M.mat')

%%

Omax=2000;
step=50;

offset= -Omax:step:Omax;

k_offset=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000]*0.5;


Data_Test_O=TMP_Zspectra_matrix;


%% B0 correction
Data_Test_CB0=B0correct_CEST(Data_Test_O,32);


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
        sig1=double(squeeze(Data_Test_CB0(ii,jj,:)));
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

        LoD_2pool(ii,jj,:)=(1./Data_Test_CB0(ii,jj,:)-1./(1-sig_simu_2pool)).*TMP_R1W_cal_matrix(ii,jj)*(1+TMP_fm_matrix(ii,jj));
        LoD_2pool_mtr(ii,jj,:) =(1-Data_Test_CB0(ii,jj,:)) - sig_simu_2pool;

        sprintf("LD----------------------- %d",ii)
    end
end