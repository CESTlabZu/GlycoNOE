function vals = matsolv_LD(pars, x, Delta)
    % matsolv computes a two-pool Lorentzian function.
    % 
    % Inputs:
    %   pars  - parameter vector [A1, b1, LW1, A2, b2, LW2]
    %   x     - vector of offset values (e.g., frequency offsets)
    %   Delta - extra parameter (not used here but kept for compatibility)
    %
    % Output:
    %   vals  - computed Lorentzian background signal
    
    A1 = pars(1);
    b1 = pars(2);
    LW1 = pars(3);
    A2 = pars(4);
    b2 = pars(5);
    LW2 = pars(6);
    
    % Compute the two Lorentzian functions and sum them:
    p1 = A1 ./ (1 + ((x - b1) ./ (0.5 * LW1)).^2);
    p2 = A2 ./ (1 + ((x - b2) ./ (0.5 * LW2)).^2);
    
    vals = p1 + p2;
end
