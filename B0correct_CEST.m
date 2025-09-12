%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BO CORRECTION BY USING LD AND SYMMETRY INFORMATION
%
% Authors: Leqi Yin, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CESTZ_B0_CORRECTED = B0correct_CEST(Zspectrum,gridsize)
% correct_CEST performs LSQ fitting and B0-correction on the noisy CEST Z-spectra.
%
% INPUT:
%   Zspectrum - 3D array of size [128 x 128 x 89] containing the noisy Z-spectra.
%
% OUTPUT:
%   CESTZ_B0_CORRECTED - 3D array [128 x 128 x 89] containing the
%                                B0-corrected Z-spectra.
%
% This function uses LSQ curve fitting with a 2-pool model to
% compute a clean version of the Z-spectrum (DL_LoD_inv) from which a B0 shift is 
% estimated (using a symmetry measure) and applied to the original noisy data.

%% Set up parameters
Omax = 1000;
step = 25;
offset = -Omax:step:Omax;
k_0 = [-2000, -1750, -1500, -1250, offset, 1250, 1500, 1750, 2000]';  % Column vector of offsets

% Parameters for the 2-pool model fitting
beta0_2pool = [0.9, 0, 280, 0.1, 0, 5000];      % initial guess
lb_2pool    = [0.02, -200, 60, 0, -800, 2000];    % lower bounds
ub_2pool    = [1, 200, 2000, 1, 800, 20000];       % upper bounds

x = k_0;  % x values for the full Z-spectrum

% Define the reduced x values for the 2-pool model (row vector of length 17)
x_2pool = [-4000; -3500; -3000; -2500; -200; -150; -100; -50; 0; 50; 100; 150; 200; 2500; 3000; 3500; 4000] / 2;
x_2pool = x_2pool';

Data = Zspectrum;  % Rename for clarity; Data is the noisy Z-spectrum

%% Preallocate the LSQ fitting output array
LoD_2pool_train = zeros(size(Data));  % same size as Data: [128 x 128 x 89]

%% Loop over each voxel for LSQ fitting
for ii = 1:gridsize
    for jj = 1:gridsize
        if any(squeeze(Data(ii,jj,:)) ~= 0)

            sig = 1 - squeeze(Data(ii,jj,:));  
            sig_2pool = zeros(1,17);
            for i = 1:4
                sig_2pool(i) = sig(i);
            end
            for i = 5:13
                sig_2pool(i) = sig(i + 36);
            end
            for i = 14:17
                sig_2pool(i) = sig(i + 72);
            end
            
            Delta = 1;  % Parameter passed to the model functions
            
            % Set LSQ curve fitting options
            options = optimset('lsqcurvefit');
            options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,...
                'MaxFunEvals',5e4*length(x),'MaxIter',2e5);
            
            % Fit the 2-pool model using lsqcurvefit
            [beta_2pool, ~, ~, ~, ~, ~, ~] = lsqcurvefit(@matsolv_2pool, beta0_2pool, ...
                x_2pool, sig_2pool, lb_2pool, ub_2pool, options, Delta);
            
            % Compute the clean Z-spectrum using the fitted parameters
            LoD_2pool_train(ii,jj,:) = matsolv_LD(beta_2pool, x, Delta);
        else
            LoD_2pool_train(ii,jj,:) = zeros(1,1,89);
        end
        
        % Display progress for each row (optional)
        disp(['LD Processing voxel (', num2str(ii), ',', num2str(jj),')']);
    end
end

% Compute the "clean" data from the LSQ fitting
DL_LoD_inv = 1 - LoD_2pool_train;

%% B0 Correction Step

% Define symmetry index pairs around the center index (45)
symPairs = [40 50; 41 49; 42 48; 43 47; 44 46];

% Define candidate shift range (adjust based on expected B0 inhomogeneity)
candidateShifts = -100:1:100;  

% Preallocate the final corrected output array
CESTZ_B0_CORRECTED = zeros(size(Data));

for ii = 1:gridsize
    for jj = 1:gridsize
        LDSpec = squeeze(DL_LoD_inv(ii,jj,:));
        Zspectrum_voxel = squeeze(Data(ii,jj,:));
        
        if any(LDSpec ~= 0)
            % Use the clean spectrum to estimate the best shift
            bestShift = findBestSymShift(LDSpec, k_0, symPairs, candidateShifts);
            
            % Apply the estimated shift to the original noisy data
            CESTZ_B0_CORRECTED(ii,jj,:) = interp1(...
                k_0 - bestShift, ...   % shifted x values
                Zspectrum_voxel, ...     % noisy data
                k_0, ...               % reference grid (keeps the same shape)
                'linear', 1);       % use 0.95 for out-of-bound values

        else
            CESTZ_B0_CORRECTED(ii,jj,:) = zeros(1,1,89);
        end
        disp(['Correcting voxel (', num2str(ii), ',', num2str(jj),')']);
    end
end


end
