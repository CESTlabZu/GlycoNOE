function [FitResult,fitParam] = plof_3(k, input, fitParam)
    input_background = input(45:85);
    k_background = k(45:85);
    FitResult.xindex = k_background;
    % FitResult.xindex = linspace(min(k_background),max(k_background),200);
    
    input_background(6:12) = [];
    k_background(6:12) = [];

    % x0_background = [0.1, 1, 0.5, -190];
    % lb=[0, 0, 0,-1000];
    % ub=[100, 100, 1000,0];
    x0_background = [10, 0.1, 0.5, -10];
    lb=[0, 0, 0,-100];
    ub=[100, 1, 1,0];
    
    options=optimset('MaxFunEvals',1e6,'TolFun',1e-6,'TolX',1e-6, 'Display',  'off' );
    [FitResult.Coefficents,resnorm]=lsqcurvefit(@CurveFunction_background,x0_background,k_background,input_background,lb,ub,options,fitParam);
    FitResult.Coefficents1=FitResult.Coefficents;
    FitResult.Background = CurveFunction_background(FitResult.Coefficents,FitResult.xindex,fitParam);


    % lb(1) = 1e-4;
    % ub(1) = 100;
    % lb(2) = 0.1;
    % ub(2) = 5;
    lb(1) = 1e-4;
    ub(1) = 0.1;
    lb(2) = 0.1;
    ub(2) = 5;
    lb(3) = min(fitParam.PeakRange);
    ub(3) = max(fitParam.PeakRange);

    lb(4:7)=FitResult.Coefficents - 0*abs(FitResult.Coefficents*0.02); % fix the background parameters
    ub(4:7)=FitResult.Coefficents + 0*abs(FitResult.Coefficents*0.02);
    x0 = [0.01, 1, 1, 0,0,0,0];
    
    [FitResult.Coefficents,resnorm]=lsqcurvefit(@CurveFunction,x0,k(45:85),input(45:85),lb,ub,options,fitParam);
    FitResult.Coefficents2=FitResult.Coefficents;
    FitResult.Curve = CurveFunction(FitResult.Coefficents,FitResult.xindex,fitParam);

    FitResult.Offset = k';
    FitResult.Saturation = input;
    FitResult.Offset_background = input_background;
    FitResult.Saturation_background = k_background;
    FitResult.Arex = ((1./FitResult.Curve)-(1./FitResult.Background)).*fitParam.R1.*(1+fitParam.fm);
    %------------ Calculate deltaZ and assign the Coefficents ---------------%
    FitResult.Rpeak = FitResult.Coefficents(1);
    FitResult.FitPeakOffset = FitResult.Coefficents(3);

end