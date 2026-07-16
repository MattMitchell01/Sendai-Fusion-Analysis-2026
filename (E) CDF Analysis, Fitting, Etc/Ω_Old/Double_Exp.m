function [AllResults,ResultsReport] = Double_Exp(CumX,CumYNorm,CumYDecayNorm,FigureHandles,...
    CurrentColor,FileNumber,ResultsReport,UsefulInfo,Options,AllResults)


NumberDataPoints = ResultsReport(FileNumber).NumVirusInCDF;
% CumYNorm = max(CumYDecayNorm) - CumYDecayNorm;

    % Set up bounds and initial guesses. [Amplitude,Tau, Tau 2] Amplitude 2 is defined by Amplitude 1
    Low_Bounds = [1e-6,1e-6,1e-6];
    Up_Bounds = [1,max(CumX)/2,max(CumX)/2];    
    Init_Guess = [0.5,(ResultsReport(FileNumber).MeanpHtoFuse)/4,ResultsReport(FileNumber).MeanpHtoFuse];

    
if strcmp(UsefulInfo.CurrentFitMethod,'Max Like')
    FitName = '---2 Exp---';
    OptimizationOptions = optimset('TolX', 1e-6,'TolFun', 1e-9,'Algorithm', 'interior-point','Display','off');

    [FitValues,LogLike, ExitFlag] = fmincon(@Double_Exp_MLE,Init_Guess,[],[],[],[],...
        Low_Bounds,Up_Bounds,[],OptimizationOptions,CumX,CumYDecayNorm,NumberDataPoints);
    LogLike = -LogLike;
elseif strcmp(UsefulInfo.CurrentFitMethod,'Least Squares')
else
    disp(' Error, wrong fit method chosen')
end

    %  Plot CDF and fit line
    if UsefulInfo.FitNumber == 1
        set(0,'CurrentFigure',FigureHandles.FitWindow)
        hold on
        plot(CumX,CumYNorm,CurrentColor.DataPoints);
    end

    Prop1 = FitValues(1);
    Tau1 = FitValues(2);
    Prop2 = 1-Prop1;
    Tau2 = FitValues(3);
    
        FitLine = Prop1 * exp(-CumX./Tau1) + Prop2 * exp(-CumX./Tau2);
                FitLine = 1 - FitLine;
%         set(0,'CurrentFigure',FigureHandles.FitWindow)
%         hold on
%         plot(CumX,FitLine,CurrentColor.FitLine,'LineWidth',2)
        
        %Make sure we plot back to zero (in the case our data starts much
        %later)
        XtoPlot = 0:1:max(CumX);
        FitLineToPlot = Prop1 * exp(-XtoPlot./Tau1) + Prop2 * exp(-XtoPlot./Tau2);
                FitLineToPlot = 1 - FitLineToPlot;
        set(0,'CurrentFigure',FigureHandles.FitWindow)
        hold on
        
        plot(XtoPlot,FitLineToPlot,CurrentColor.FitLine,'LineWidth',2)
        
        %R squared calculation
        MeanData = mean(CumYNorm);
        SST = sum((CumYNorm-MeanData).^2);
        SSE = sum((CumYNorm-FitLine).^2);
        R2 = 1-SSE/SST;
        Residuals = CumYNorm-FitLine;

        %Plot residuals
        set(0,'CurrentFigure',FigureHandles.ResidualsWindow)
        hold on
        plot(CumX,Residuals,CurrentColor.ResidualPoints);

        %AIC calculation
            NumberofParameters = length(FitValues);
            NumberPointsinFit = length(CumX);
            AIC = 2*NumberofParameters + 2*LogLike;
%                 [AIC, BIC] = aicbic(LogLike,NumberofParameters,NumberPointsinFit);
%                 AICC = AIC - (2*NumberofParameters*(NumberofParameters+1)/(NumberPointsinFit-NumberofParameters-1));


                ParameterMatrix = cell(2,length(FitValues));
                    ParameterMatrix{1,1} ='Prop1';
                    ParameterMatrix{1,2} ='Tau1';
                    ParameterMatrix{1,3} ='Prop2';
                    ParameterMatrix{1,4} ='Tau2';
                    ParameterMatrix{1,5} ='R2';
                    
                    ParameterMatrix{2,1} =Prop1;
                    ParameterMatrix{2,2} =Tau1;
                    ParameterMatrix{2,3} =1-Prop1;
                    ParameterMatrix{2,4} =Tau2;
                    ParameterMatrix{2,5} =R2;
                    
                FitStatsMatrix = cell(1,3);
                    FitStatsMatrix{1,1} =AIC;
                    FitStatsMatrix{1,2} =LogLike;
                    FitStatsMatrix{1,3} =R2;
                    
                FitNumber = UsefulInfo.FitNumber;
                    AllResults(FileNumber).FitResults(FitNumber).ParameterMatrix = ParameterMatrix;
                    AllResults(FileNumber).FitResults(FitNumber).FitStatsMatrix = FitStatsMatrix;
                    AllResults(FileNumber).FitResults(FitNumber).FitType=FitName;
                    AllResults(FileNumber).FitResults(FitNumber).Prop1=FitValues(1);
                    AllResults(FileNumber).FitResults(FitNumber).Tau1=FitValues(2);
                    AllResults(FileNumber).FitResults(FitNumber).Prop2=Prop2;
                    AllResults(FileNumber).FitResults(FitNumber).Tau2=Tau2;
                    AllResults(FileNumber).FitResults(FitNumber).LogLike= LogLike;
                    AllResults(FileNumber).FitResults(FitNumber).AIC= AIC;
                    AllResults(FileNumber).FitResults(FitNumber).R2=R2;
%                     
%         if UsefulInfo.FitNumber == 1
%             ResultsReport(FileNumber).FitType= FitName;
%             ResultsReport(FileNumber).Amp1_1=FitValues(1);
%             ResultsReport(FileNumber).Tau1_1=FitValues(2);
%             ResultsReport(FileNumber).Amp2_1=Prop2;
%             ResultsReport(FileNumber).Tau2_1=Tau2;
%             ResultsReport(FileNumber).LogLike1 = LogLike;
%             ResultsReport(FileNumber).AIC1 = AIC;
%             ResultsReport(FileNumber).R2=R2;
%         elseif UsefulInfo.FitNumber == 2
%             ResultsReport(FileNumber).FitType2=FitName;
%             ResultsReport(FileNumber).Amp1_2=FitValues(1);
%             ResultsReport(FileNumber).Tau1_2=FitValues(2);
%             ResultsReport(FileNumber).Amp2_2=Prop2;
%             ResultsReport(FileNumber).Tau2_2=Tau2;
%             ResultsReport(FileNumber).LogLike2= LogLike;
%             ResultsReport(FileNumber).AIC2 = AIC;
%             ResultsReport(FileNumber).R2_2=R2;
%         elseif UsefulInfo.FitNumber == 3
%             ResultsReport(FileNumber).FitType3=FitName;
%             ResultsReport(FileNumber).Amp1_3=FitValues(1);
%             ResultsReport(FileNumber).Tau1_3=FitValues(2);
%             ResultsReport(FileNumber).Amp2_3=Prop2;
%             ResultsReport(FileNumber).Tau2_3=Tau2;
%             ResultsReport(FileNumber).LogLike3= LogLike;
%             ResultsReport(FileNumber).AIC3 = AIC;
%             ResultsReport(FileNumber).R2_3=R2;
%         elseif UsefulInfo.FitNumber == 4
%             ResultsReport(FileNumber).FitType4=FitName;
%             ResultsReport(FileNumber).Amp1_4=FitValues(1);
%             ResultsReport(FileNumber).Tau1_4=FitValues(2);
%             ResultsReport(FileNumber).Amp2_4=Prop2;
%             ResultsReport(FileNumber).Tau2_4=Tau2;
%             ResultsReport(FileNumber).LogLike4= LogLike;
%             ResultsReport(FileNumber).AIC4 = AIC;
%             ResultsReport(FileNumber).R2_4=R2;
%         else
%             disp(' Error: you chose too many fits to compare')
%         end
%         
end


function LogLike = Double_Exp_MLE(p,CumX,CumY,NumberDataPoints)
    x=CumX;
    n=NumberDataPoints;
    y=n*CumY;
    
    %Model
    Prob = p(1)*exp(-x./p(2))+(1 - p(1))*exp(-x./p(3));
        Prob = Prob + (Prob == 0).*1e-5 - (Prob==1).*1e-5; %Ensure 0<Prob<1
    
    %Log Likelihood
    LogLike = (-1)*(y.*log(Prob)+(n-y).*log(1-Prob));
    LogLike = sum(LogLike);
end