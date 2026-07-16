function [AllDataToAnalyze,CDFReports,FitLine] = Fit_Type_Sorter(AllDataToAnalyze,CDFReports, Options, FitType, FileNumber)
    
    if strcmp(FitType,'1 Exp')
        [CDFReports,FitLine] = Model_Single_Exp(AllDataToAnalyze,CDFReports, Options, FileNumber);
    elseif strcmp(FitType,'Gamma')
        [CDFReports,FitLine] = Model_Gamma(AllDataToAnalyze,CDFReports, Options, FileNumber);
    else
        disp("!!!!ERROR: Set up your fit type correctly please.")
        StopProgramNow
    end
    

end