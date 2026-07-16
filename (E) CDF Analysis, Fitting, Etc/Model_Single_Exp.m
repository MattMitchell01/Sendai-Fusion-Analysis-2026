function [CDFReports,FitLine] = Model_Single_Exp(AllDataToAnalyze,CDFReports, Options, FileNumber)

    CurrCDFData = AllDataToAnalyze(FileNumber).CDFData;
    CumXToFit = CurrCDFData.CumX;
    CumYToFit = CurrCDFData.CumYNorm;

    if strcmp(Options.AccountForDeadTime,'y') 
        DeadTime = CurrCDFData.DeadTime;
    else
        DeadTime = NaN;
    end

    % Set up bounds and initial guesses. [Tau]
        Init_Guess = mean(CurrCDFData.SortedWaitTimeList);
        Low_Bounds = 0;
        Up_Bounds = Init_Guess*10;
    
        OptimizationOptions = optimset('TolX', 1e-6,'TolFun', 1e-9,'Algorithm', 'interior-point','Display','off');
    
    % Run fit
        [FitValues,SSE] = fmincon(@Single_Exp_LSQ_Model,Init_Guess,[],[],[],[],...
            Low_Bounds,Up_Bounds,[],OptimizationOptions,CumXToFit,CumYToFit,DeadTime);

    % Generate fit line to plot
        FitXValues = CumXToFit;
        tau = FitValues(1);

        if isnan(DeadTime)
            FitYValues = expcdf(CumXToFit,tau);
            ExtrapXZeroVal = 0;
        else
            DeadYVal = expcdf(DeadTime,tau);
            FitYValues = (expcdf(CumXToFit,tau)-DeadYVal)./(1-DeadYVal);
            ExtrapXZeroVal = (expcdf(0,tau)-DeadYVal)./(1-DeadYVal);
        end

        FitLine.FitXValues = FitXValues;
        FitLine.FitYValues = FitYValues;
 
    % AIC calc 
    % (based on https://en.wikipedia.org/wiki/Akaike_information_criterion#Comparison_with_least_squares)
        NumFitParams = 1; %for single exp
        NumDataPts = length(CurrCDFData.SortedWaitTimeList);
        MLEEstimator = SSE/NumDataPts;
        AIC = 2*NumFitParams + NumDataPts*log(MLEEstimator);

    % Add to results report
        CDFReports(FileNumber).Exp_FitName = 'Fit: 1 Exp';
        CDFReports(FileNumber).ExpFit_Tau = tau;
        CDFReports(FileNumber).ExpFit_AIC = AIC;
        CDFReports(FileNumber).ExpFit_DeadTimeAccountedFor = Options.AccountForDeadTime;
        if strcmp(Options.AccountForDeadTime,'y') 
            CDFReports(FileNumber).ExpFit_YValAtXZero = ExtrapXZeroVal;
        end

end

function [SSE] = Single_Exp_LSQ_Model(p,CumXToFit,CumYToFit,DeadTime)

    tau =p(1);

    if isnan(DeadTime)
        FitYValues = expcdf(CumXToFit,tau);
    else
        DeadYVal = expcdf(DeadTime,tau);
        FitYValues = (expcdf(CumXToFit,tau)-DeadYVal)./(1-DeadYVal);
    end
    
    SSE = sum((CumYToFit - FitYValues).^2);

end

