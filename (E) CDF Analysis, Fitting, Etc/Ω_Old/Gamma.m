function [AllResults,ResultsReport] = Gamma(CumX,CumYNorm,FigureHandles,CurrentColor,FileNumber,ResultsReport,UsefulInfo,Options,AllResults)

NumberDataPoints = ResultsReport(FileNumber).NumVirusInCDF;

if strcmp(UsefulInfo.CurrentFitMethod,'Max Like')
    OptimizationOptions = optimset('TolX', 1e-6,'TolFun', 1e-9,'Algorithm', 'interior-point');
    
    
    % Set up bounds and initial guesses. [N,Tau]
    [FitName,LowBounds,UpBounds,InitialGuess] = Setup_Fit_Parameters_Gamma(Options,CumX,ResultsReport,FileNumber,UsefulInfo);

    [FitValues,LogLike, ExitFlag] = fmincon(@Gamma_MLE_Model,InitialGuess,[],[],[],[],...
        LowBounds,UpBounds,[],OptimizationOptions,CumX,CumYNorm,NumberDataPoints);

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

                NumberSteps = FitValues(1);
                Tau = FitValues(2);
                FitLine = gamcdf(CumX,NumberSteps,Tau);
%                 set(0,'CurrentFigure',FigureHandles.FitWindow)
%                 hold on
%                 plot(CumX,FitLine,CurrentColor.FitLine,'LineWidth',2)

                %Make sure we plot back to zero (in the case our data starts much
                %later)
                XtoPlot = 0:1:max(CumX);
                FitLineToPlot = gamcdf(XtoPlot,NumberSteps,Tau);
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
                if  strcmp(Options.FixGammaN,'y')
                    NumberofParameters = length(FitValues) - 1;
                else
                    NumberofParameters = length(FitValues);
                end
                
                NumberPointsinFit = length(CumX);
                AIC = 2*NumberofParameters + 2*LogLike;
%                 [AIC, BIC] = aicbic(LogLike,NumberofParameters,NumberPointsinFit);
%                 AICC = AIC - (2*NumberofParameters*(NumberofParameters+1)/(NumberPointsinFit-NumberofParameters-1));

                FitNumber = UsefulInfo.FitNumber;
                    AllResults(FileNumber).FitResults(FitNumber).FitType=FitName;
                    AllResults(FileNumber).FitResults(FitNumber).NumberSteps=FitValues(1);
                    AllResults(FileNumber).FitResults(FitNumber).Tau=FitValues(2);
                    AllResults(FileNumber).FitResults(FitNumber).LogLike= LogLike;
                    AllResults(FileNumber).FitResults(FitNumber).AIC= AIC;
                    AllResults(FileNumber).FitResults(FitNumber).R2=R2;
                    
                
                if UsefulInfo.FitNumber == 1
                    ResultsReport(FileNumber).FitType1=FitName;
                    ResultsReport(FileNumber).NumberSteps1=FitValues(1);
                    ResultsReport(FileNumber).Tau1 =FitValues(2);
                    ResultsReport(FileNumber).LogLike1 = LogLike;
                    ResultsReport(FileNumber).AIC1 = AIC;
                    ResultsReport(FileNumber).R2=R2;
                elseif UsefulInfo.FitNumber == 2
                    ResultsReport(FileNumber).FitType2=FitName;
                    ResultsReport(FileNumber).NumberSteps2=FitValues(1);
                    ResultsReport(FileNumber).Tau2=FitValues(2);
                    ResultsReport(FileNumber).LogLike2= LogLike;
                    ResultsReport(FileNumber).AIC2 = AIC;
                    ResultsReport(FileNumber).R2_2=R2;
                elseif UsefulInfo.FitNumber == 3
                    ResultsReport(FileNumber).FitType3=FitName;
                    ResultsReport(FileNumber).NumberSteps3=FitValues(1);
                    ResultsReport(FileNumber).Tau3=FitValues(2);
                    ResultsReport(FileNumber).LogLike3= LogLike;
                    ResultsReport(FileNumber).AIC3 = AIC;
                    ResultsReport(FileNumber).R2_3=R2;
                elseif UsefulInfo.FitNumber == 4
                    ResultsReport(FileNumber).FitType4=FitName;
                    ResultsReport(FileNumber).NumberSteps4=FitValues(1);
                    ResultsReport(FileNumber).Tau4=FitValues(2);
                    ResultsReport(FileNumber).LogLike4= LogLike;
                    ResultsReport(FileNumber).AIC4 = AIC;
                    ResultsReport(FileNumber).R2_4=R2;
                else
                    disp(' Error: you chose too many fits to compare')
                end
                
                AllResults(FileNumber).ResultsReport = ResultsReport(FileNumber);
                
%             FitNumber = UsefulInfo.FitNumber;
%             
%                     ResultsReport(FileNumber).GammaFitResults(FitNumber).FitType=FitName;
%                     ResultsReport(FileNumber).GammaFitResults(FitNumber).NumberSteps=FitValues(1);
%                     ResultsReport(FileNumber).GammaFitResults(FitNumber).Tau=FitValues(2);
%                     ResultsReport(FileNumber).GammaFitResults(FitNumber).LogLike= LogLike;
%                     ResultsReport(FileNumber).GammaFitResults(FitNumber).AIC= AIC;
%                     ResultsReport(FileNumber).GammaFitResults(FitNumber).R2=R2;
end

function LogLike = Gamma_MLE_Model(p,CumX,CumY,NumberDataPoints)
    x=CumX;
    n=NumberDataPoints;
    y=n*CumY;
    NumberSteps = p(1);
    Tau = p(2);
    
    %Model
    Prob = gamcdf(x,NumberSteps,Tau);
        %Prob = Prob + (Prob == 0).*1e-5 - (Prob==1).*1e-5; %Ensure 0<Prob<1
    
    %Log Likelihood
    LogLike = (-1)*(y.*log(Prob)+(n-y).*log(1-Prob));
    LogLike = sum(LogLike);
end

function [FitName,LowBounds,UpBounds,InitialGuess] = Setup_Fit_Parameters_Gamma(Options,CumX,ResultsReport,FileNumber,UsefulInfo)

FitNumber = UsefulInfo.FitNumber;

% Set up bounds and initial guesses. [N,Tau]
    
    if strcmp(Options.FixGammaN,'y')
        FixedValue = Options.FixedNValues(FitNumber);
        
        LowBounds = [FixedValue - 0.01,1];
        UpBounds = [FixedValue + 0.01,max(CumX)];
        InitialGuess = [5,ResultsReport(FileNumber).MeanpHtoFuse/2];
        
        FitName =  strcat('--Gamma MLE, N=', num2str(FixedValue),'--');
    else 
        FitName = '----Gamma MLE----';
        LowBounds = [1,1];
        UpBounds = [10,max(CumX)];
        InitialGuess = [5,ResultsReport(FileNumber).MeanpHtoFuse/2];
    end
end