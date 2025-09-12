function vals = matsolv_Voigt(pars, x, Delta)

A1=pars(1);
b1=pars(2);
LW1=pars(3);
A2=pars(4);
b2=pars(5);
LW2=pars(6);
A3=pars(7);
b3=pars(8);
LW3=pars(9);
A4=pars(10);
b4=pars(11);
LW4=pars(12);
A5=pars(13);
b5=pars(14);
LW5=pars(15);
A6=pars(16);
b6=pars(17);
LW6=pars(18);


off=x;



% p=A1*1./(1+((off-b1)./(0.5.*LW1)).^2)+ A2*1./(1+((off-b2)./(0.5.*LW2)).^2)+A3*1./(1+((off-b3)./(0.5.*LW3)).^2)+A4*1./(1+((off-b4)./(0.5.*LW4)).^2)+A5*1./(1+((off-b5)./(0.5.*LW5)).^2)+A6*1./(1+((off-b6)./(0.5.*LW6)).^2);
% Convert Gaussian FWHM (Delta) to standard deviation
    sigma = Delta / (2 * sqrt(2 * log(2))); % Standard deviation of Gaussian

    % Precompute Lorentzian half-widths
    gamma1 = LW1 / 2; gamma2 = LW2 / 2; gamma3 = LW3 / 2;
    gamma4 = LW4 / 2; gamma5 = LW5 / 2; gamma6 = LW6 / 2;

    % Precompute Lorentzian weights
    eta1 = gamma1 / (gamma1 + sigma * sqrt(2));
    eta2 = gamma2 / (gamma2 + sigma * sqrt(2));
    eta3 = gamma3 / (gamma3 + sigma * sqrt(2));
    eta4 = gamma4 / (gamma4 + sigma * sqrt(2));
    eta5 = gamma5 / (gamma5 + sigma * sqrt(2));
    eta6 = gamma6 / (gamma6 + sigma * sqrt(2));

    % Compute Voigt profile contributions for each peak
    V1 = A1 * (eta1 ./ (1 + ((off - b1) ./ (0.5 * LW1)).^2) + ...
        (1 - eta1) * exp(-((off - b1).^2) / (2 * sigma^2)));
    V2 = A2 * (eta2 ./ (1 + ((off - b2) ./ (0.5 * LW2)).^2) + ...
        (1 - eta2) * exp(-((off - b2).^2) / (2 * sigma^2)));
    V3 = A3 * (eta3 ./ (1 + ((off - b3) ./ (0.5 * LW3)).^2) + ...
        (1 - eta3) * exp(-((off - b3).^2) / (2 * sigma^2)));
    V4 = A4 * (eta4 ./ (1 + ((off - b4) ./ (0.5 * LW4)).^2) + ...
        (1 - eta4) * exp(-((off - b4).^2) / (2 * sigma^2)));
    V5 = A5 * (eta5 ./ (1 + ((off - b5) ./ (0.5 * LW5)).^2) + ...
        (1 - eta5) * exp(-((off - b5).^2) / (2 * sigma^2)));
    V6 = A6 * (eta6 ./ (1 + ((off - b6) ./ (0.5 * LW6)).^2) + ...
        (1 - eta6) * exp(-((off - b6).^2) / (2 * sigma^2)));

    % Sum contributions from all peaks
    vals = V1 + V2 + V3 + V4 + V5 + V6;

end