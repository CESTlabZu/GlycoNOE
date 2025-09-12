function [y, Rbak] = VoigtFunction_background(params, xdata, FitParam)
    % Parameters:
    % params: Optimization parameters
    %   params(1): Amplitude of the Lorentzian component
    %   params(2): Width (HWHM) of the Lorentzian component
    %   params(3): Constant background offset
    %   params(4): Linear background slope
    % xdata: Independent variable (frequency offsets)
    % FitParam: Additional fixed parameters (e.g., saturation power, offsets)

    % Extract the peak offset from FitParam
    peakoffset = FitParam.PeakOffset;

    % Lorentzian component for background
    Lorentzian = params(1) * params(2)^2 ./ (params(2)^2 + 4 * xdata.^2);

    % Gaussian component for background
    Gaussian = params(1) * exp(-((xdata).^2) / (2 * (params(2)^2)));

    % Combine Gaussian and Lorentzian with a weighting factor
    eta = 0.5;  % Mixing factor (0 = pure Gaussian, 1 = pure Lorentzian)
    Rbak = eta * Lorentzian + (1 - eta) * Gaussian + ...
           params(3) + params(4) * 0.001 * (xdata - peakoffset);

    % Saturation power in Hz
    satHz = FitParam.satpwr * 42.58;

    % Magnetic field effect
    cos2thet = 1 - satHz^2 ./ (satHz^2 + (FitParam.Magfield * xdata).^2);

    % Total relaxation rate
    Rall = Rbak;

    % Final background signal model
    % y = (1 - cos2thet * FitParam.R1 ./ Rall) .* exp(-Rall * FitParam.tsat) + ...
        % cos2thet * FitParam.R1 ./ Rall;

    y = cos2thet .* FitParam.R1 ./ Rall;

end
