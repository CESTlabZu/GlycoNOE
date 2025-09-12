function y = CurveFunction(x, xdata, FitParam)
    % Parameters:
    % x(1): Amplitude of the Lorentzian component for PCr
    % x(2): Width (HWHM) of the Lorentzian component for PCr
    % x(3): Offset of the PCr peak
    % x(4): Amplitude of the Lorentzian component for the background
    % x(5): Width (HWHM) of the Lorentzian component for the background
    % x(6): Constant background offset
    % x(7): Linear background slope

    % PCr resonance (Lorentzian)
    RPCr = x(1) * x(2)^2 ./ (x(2)^2 + 4 * (xdata - x(3)).^2);

    % Background signal using pseudo-Voigt approach
    % Gaussian component
    Gaussian = x(4) * exp(-((xdata).^2) / (2 * x(5)^2));

    % Lorentzian component
    Lorentzian = x(4) * x(5)^2 ./ (x(5)^2 + 4 * xdata.^2);

    % Mixing factor (eta)
    eta = 0.5; % Adjust or make dynamic if needed

    % Combined background signal
    Rbak = eta * Lorentzian + (1 - eta) * Gaussian + ...
           x(6) + x(7) * 0.001 * (xdata - FitParam.PeakOffset);

    % Saturation power in Hz
    satHz = FitParam.satpwr * 42.58;

    % Magnetic field effect
    cos2thet = 1 - satHz^2 ./ (satHz^2 + (FitParam.Magfield * xdata).^2);

    % Total relaxation rate
    Rall = Rbak + RPCr;

    % Final signal model
    % y = (1 - cos2thet * FitParam.R1 ./ Rall) .* exp(-Rall * FitParam.tsat) + ...
    %     cos2thet * FitParam.R1 ./ Rall;
    y = cos2thet .* FitParam.R1 ./ Rall;
end
