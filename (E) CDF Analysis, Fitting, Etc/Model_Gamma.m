function [CDFReports,FitLine] = Model_Gamma(AllDataToAnalyze,CDFReports, Options, FileNumber)

    CurrCDFData = AllDataToAnalyze(FileNumber).CDFData;
    CumXToFit = CurrCDFData.CumX;
    CumYToFit = CurrCDFData.CumYNorm;

    if strcmp(Options.AccountForDeadTime,'y') 
        DeadTime = CurrCDFData.DeadTime;
    else
        DeadTime = NaN;
    end

    % Set up bounds and initial guesses. Order = [NGammaSteps,Tau]
        MeanWaitTime = mean(CurrCDFData.SortedWaitTimeList);
        Init_Guess = [2,MeanWaitTime];
        Low_Bounds = [0,0];
        Up_Bounds = [10,MeanWaitTime*10];
    
        OptimizationOptions = optimset('TolX', 1e-6,'TolFun', 1e-9,'Algorithm', 'interior-point','Display','off');
    
    % Run fit
        [FitValues, SSE] = fmincon(@Gamma_LSQ_Model,Init_Guess,[],[],[],[],...
            Low_Bounds,Up_Bounds,[],OptimizationOptions,CumXToFit,CumYToFit,DeadTime);

    % Generate fit line to plot
        FitXValues = CumXToFit;
        NGammaSteps = FitValues(1);
        tau =FitValues(2);
        
        if isnan(DeadTime)
            FitYValues = gamcdf(FitXValues,NGammaSteps,tau);
            ExtrapXZeroVal = 0;
        else
            DeadYVal = gamcdf(DeadTime,NGammaSteps,tau);
            FitYValues = (gamcdf(FitXValues,NGammaSteps,tau)-DeadYVal)./(1-DeadYVal);
            ExtrapXZeroVal = (gamcdf(0,tau)-DeadYVal)./(1-DeadYVal);
        end


        FitLine.FitXValues = FitXValues;
        FitLine.FitYValues = FitYValues;

    % AIC calc 
    % (based on https://en.wikipedia.org/wiki/Akaike_information_criterion#Comparison_with_least_squares)
        NumFitParams = 2; %for gamma
        NumDataPts = length(CurrCDFData.SortedWaitTimeList);
        MLEEstimator = SSE/NumDataPts;
        AIC = 2*NumFitParams + NumDataPts*log(MLEEstimator);

    % Add to results report
        CDFReports(FileNumber).Gamma_FitName = 'Fit: Gamma';
        CDFReports(FileNumber).GammaFit_NSteps = NGammaSteps;
        CDFReports(FileNumber).GammaFit_Tau = tau;
        CDFReports(FileNumber).GammaFit_AIC = AIC;
        CDFReports(FileNumber).GammaFit_DeadTimeAccountedFor = Options.AccountForDeadTime;
        if strcmp(Options.AccountForDeadTime,'y') 
            CDFReports(FileNumber).GammaFit_YValAtXZero = ExtrapXZeroVal;
        end

end

function [SSE] = Gamma_LSQ_Model(p,CumXToFit,CumYToFit,DeadTime)

    NGammaSteps = p(1);
    tau =p(2);

    if isnan(DeadTime)
        FitYValues = gamcdf(CumXToFit,NGammaSteps,tau);
    else
        DeadYVal = gamcdf(DeadTime,NGammaSteps,tau);
        FitYValues = (gamcdf(CumXToFit,NGammaSteps,tau)-DeadYVal)./(1-DeadYVal);
    end
    
    SSE = sum((CumYToFit - FitYValues).^2);

end

