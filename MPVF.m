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


%% multiple-pool Voigt fit.

for i=1:32
     for j=1:32
        
        sig=(1-double(squeeze(Data_Test_CB0(i,j,:))));
        
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
        beta_matrix_MPVF(i,j,:)=beta;
        sig_simur_GCG=matsolv_Voigt(beta_GCG,x,Delta);
        beta_GCG(10)=0;
        sig_simur_ref_GCG=matsolv_Voigt(beta_GCG,x,Delta);



        MPVF_AREX_GCG(i,j,:)=(1./(1-sig_simur_GCG)-1./(1-sig_simur_ref_GCG))*R1W_AREX*(1+fm_AREX);%use this one!!!


       
        sprintf("MPVF----------------------- %d",i)
            
     end
end