function [AllResults,ResultsReport] = Double_Gamma(CumX,CumYNorm,FigureHandles,CurrentColor,FileNumber,ResultsReport,UsefulInfo,Options,AllResults)

NumberDataPoints = ResultsReport(FileNumber).NumVirusInCDF;

if strcmp(UsefulInfo.CurrentFitMethod,'Max Like')
    OptimizationOptions = optimset('TolX', 1e-6,'TolFun', 1e-9,'Algorithm', 'interior-point');
    
    

    [FitName,LowBounds,UpBounds,InitialGuess] = Setup_Fit_Parameters_Double_Gamma(Options,CumX,ResultsReport,FileNumber,UsefulInfo);

    [FitValues,LogLike, ExitFlag] = fmincon(@Double_Gamma_MLE_Model,InitialGuess,[],[],[],[],...
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

                Prop1 = FitValues(1);
                NumberSteps1 = FitValues(2);
                Tau1= FitValues(3);
                NumberSteps2= FitValues(4);
                Tau2= FitValues(5);
%                 FitLine = zeros(size(CumX));
%                 IndextoUse = CumX > FitValues(3);
%                 CumXLag = CumX -LagTime;
                
%                 FitLine(IndextoUse) = gamcdf(CumX(IndextoUse) -LagTime,NumberSteps,Tau);
                FitLine = Prop1 * gamcdf(CumX,NumberSteps1,Tau1) + (1-Prop1) * gamcdf(CumX,NumberSteps2,Tau2);
%                 set(0,'CurrentFigure',FigureHandles.FitWindow)
%                 hold on
%                 plot(CumX,FitLine,CurrentColor.FitLine,'LineWidth',2)

                %Make sure we plot back to zero (in the case our data starts much
                %later)
                XtoPlot = 0:1:max(CumX);
                FitLineToPlot = Prop1 * gamcdf(XtoPlot,NumberSteps1,Tau1) + (1-Prop1) * gamcdf(XtoPlot,NumberSteps2,Tau2);
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

                ParameterMatrix = cell(2,length(FitValues));
                    ParameterMatrix{1,1} ='Prop1';
                    ParameterMatrix{1,2} ='N1';
                    ParameterMatrix{1,3} ='Tau1';
                    ParameterMatrix{1,4} ='Prop2';
                    ParameterMatrix{1,5} ='N2';
                    ParameterMatrix{1,6} ='Tau2';
                    ParameterMatrix{1,7} ='R2';
                    
                    ParameterMatrix{2,1} =Prop1;
                    ParameterMatrix{2,2} =NumberSteps1;
                    ParameterMatrix{2,3} =Tau1;
                    ParameterMatrix{2,4} =1-Prop1;
                    ParameterMatrix{2,5} =NumberSteps2;
                    ParameterMatrix{2,6} =Tau2;
                    ParameterMatrix{2,7} =R2;
                    
                FitStatsMatrix = cell(1,3);
                    FitStatsMatrix{1,1} =AIC;
                    FitStatsMatrix{1,2} =LogLike;
                    FitStatsMatrix{1,3} =R2;
                    
                FitNumber = UsefulInfo.FitNumber;
                AllResults(FileNumber).FitResults(FitNumber).ParameterMatrix = ParameterMatrix;
                AllResults(FileNumber).FitResults(FitNumber).FitStatsMatrix = FitStatsMatrix;
                AllResults(FileNumber).FitResults(FitNumber).FitType=FitName;
                AllResults(FileNumber).FitResults(FitNumber).Prop1=FitValues(1);
                AllResults(FileNumber).FitResults(FitNumber).NumberSteps1=FitValues(2);
                AllResults(FileNumber).FitResults(FitNumber).Tau1=FitValues(3);
                AllResults(FileNumber).FitResults(FitNumber).NumberSteps2=FitValues(4);
                AllResults(FileNumber).FitResults(FitNumber).Tau2=FitValues(5);
                AllResults(FileNumber).FitResults(FitNumber).LogLike= LogLike;
                AllResults(FileNumber).FitResults(FitNumber).AIC= AIC;
                AllResults(FileNumber).FitResults(FitNumber).R2=R2;
                
%                 if UsefulInfo.FitNumber == 1
%                     ResultsReport(FileNumber).FitType1=FitName;
%                     ResultsReport(FileNumber).Prop1_1=FitValues(1);
%                     ResultsReport(FileNumber).NumberSteps1_1=FitValues(2);
%                     ResultsReport(FileNumber).Tau1_1=FitValues(3);
%                     ResultsReport(FileNumber).NumberSteps2_1=FitValues(4);
%                     ResultsReport(FileNumber).Tau2_1=FitValues(5);
%                     ResultsReport(FileNumber).LogLike1 = LogLike;
%                     ResultsReport(FileNumber).AIC1 = AIC;
%                     ResultsReport(FileNumber).R2=R2;
%                 elseif UsefulInfo.FitNumber == 2
%                     ResultsReport(FileNumber).FitType2=FitName;
%                     ResultsReport(FileNumber).Prop1_2=FitValues(1);
%                     ResultsReport(FileNumber).NumberSteps1_2=FitValues(2);
%                     ResultsReport(FileNumber).Tau1_2=FitValues(3);
%                     ResultsReport(FileNumber).NumberSteps2_2=FitValues(4);
%                     ResultsReport(FileNumber).Tau2_2=FitValues(5);
%                     ResultsReport(FileNumber).LogLike2= LogLike;
%                     ResultsReport(FileNumber).AIC2 = AIC;
%                     ResultsReport(FileNumber).R2_2=R2;
%                 elseif UsefulInfo.FitNumber == 3
%                     ResultsReport(FileNumber).FitType3=FitName;
%                     ResultsReport(FileNumber).NumberSteps3=FitValues(1);
%                     ResultsReport(FileNumber).Tau3=FitValues(2);
%                     ResultsReport(FileNumber).LagTime3 =FitValues(3);
%                     ResultsReport(FileNumber).LogLike3= LogLike;
%                     ResultsReport(FileNumber).AIC3 = AIC;
%                     ResultsReport(FileNumber).R2_3=R2;
%                 elseif UsefulInfo.FitNumber == 4
%                     ResultsReport(FileNumber).FitType4=FitName;
%                     ResultsReport(FileNumber).NumberSteps4=FitValues(1);
%                     ResultsReport(FileNumber).Tau4=FitValues(2);
%                     ResultsReport(FileNumber).LagTime4 =FitValues(3);
%                     ResultsReport(FileNumber).LogLike4= LogLike;
%                     ResultsReport(FileNumber).AIC4 = AIC;
%                     ResultsReport(FileNumber).R2_4=R2;
%                 else
%                     disp(' Error: you chose too many fits to compare')
%                 end
%                 
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

function LogLike = Double_Gamma_MLE_Model(p,CumX,CumY,NumberDataPoints)
    x=CumX;
    n=NumberDataPoints;
    y=n*CumY;
    Prop1 = p(1);
    NumberSteps1 = p(2);
    Tau1= p(3);
    NumberSteps2= p(4);
    Tau2= p(5);
    
%     x =x -LagTime;
%     Prob = zeros(size(x))+1e-5;
%     IndextoUse = x>LagTime;
    %Model
%     Prob(IndextoUse) = gamcdf(x(IndextoUse)-LagTime,NumberSteps1,Tau1);
    Prob = Prop1 * gamcdf(x,NumberSteps1,Tau1) + (1-Prop1) * gamcdf(x,NumberSteps2,Tau2);
        Prob = Prob + (Prob == 0).*1e-5 - (Prob==1).*1e-5; %Ensure 0<Prob<1
    
    %Log Likelihood
    LogLike = (-1)*(y.*log(Prob)+(n-y).*log(1-Prob));
    LogLike = sum(LogLike);
end

function [FitName,LowBounds,UpBounds,InitialGuess] = Setup_Fit_Parameters_Double_Gamma(Options,CumX,ResultsReport,FileNumber,UsefulInfo)

FitNumber = UsefulInfo.FitNumber;

% Set up bounds and initial guesses. [Prop1,N1,Tau1,N2,Tau2] Prop2 is defined by Prop1
    
    if strcmp(Options.FixGammaN,'y')
        FixedValue = Options.FixedNValues(FitNumber);
        
        LowBounds = [0,FixedValue - 0.01,1,FixedValue - 0.01,1];
        UpBounds = [1,FixedValue + 0.01,max(CumX),FixedValue + 0.01,max(CumX)];
        InitialGuess = [ 0.5, FixedValue,ResultsReport(FileNumber).MeanpHtoFuse/2,FixedValue,ResultsReport(FileNumber).MeanpHtoFuse];
        
        FitName =  strcat('2 Gamma, N=', num2str(FixedValue));
    else 
        FitName = '2 Gamma';
        LowBounds = [0,0.1,1,0.1,1];
        UpBounds = [1,10,max(CumX),10,max(CumX)];
        InitialGuess = [0.5,2,ResultsReport(FileNumber).MeanpHtoFuse/2,5,ResultsReport(FileNumber).MeanpHtoFuse];
    end
end