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


%%

%Hybrid(LD+MPLF)

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