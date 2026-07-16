function [FitParamMatrix,GroupedResultsMatrix,NumberMatrices] = Display_Results_Matrix(Results,UsefulInfo,Options)

SkipResultsMatrix = 'n';
MinNameLength = 20;

% ResultsOption = 'Randomness';
if UsefulInfo.NumberFitsToPerform == 1 && strcmp(Options.FitTypes,'Gamma W Lag')
    ResultsOption = {'Randomness','Single Fit'};
    % ResultsOption = 'Single Fit';
elseif strcmp(Options.FixGammaN,'y') 
    ResultsOption = {'Randomness','Constrained Best Fit'};
elseif UsefulInfo.NumberFitsToPerform > 1
    ResultsOption = {'Randomness','Multiple Fit, AIC'};
else
%     SkipResultsMatrix = 'y';
    ResultsOption = 'Randomness';
end
% ResultsOption = 'Fit Stats';
% ResultsOption = 'Randomness';


if SkipResultsMatrix == 'y'
    GroupedResultsMatrix = [];
    NumberMatrices = [];
    FitParamMatrix(1).ResultsMatrix = [];
else

    if iscell(ResultsOption)
        NumberMatrices = length(ResultsOption);
    else
        NumberMatrices = 1;
    end

    for w = 1:NumberMatrices
        if iscell(ResultsOption)
            CurrentOption = ResultsOption{1,w};
        else 
            CurrentOption = ResultsOption;
        end

        [GroupedResultsMatrix(w).ResultsMatrix] = Generate_Results_Matrix(Results,UsefulInfo,CurrentOption,MinNameLength);
    end
    
    [FitParamMatrix] = Generate_Fit_Parameters_Matrix(Results,UsefulInfo,ResultsOption,MinNameLength);
end
 
end

function [FitParamMatrix] = Generate_Fit_Parameters_Matrix(Results,UsefulInfo,ResultsOption,MinNameLength)
    NumberFiles = length(Results);

    NumberFits = UsefulInfo.NumberFitsToPerform;
    
    
    for h = 1:NumberFits
        NumberParameters = size(Results(1).FitResults(h).ParameterMatrix,2)+2;
        ResultsMatrix = cell(NumberFiles +1,NumberParameters + 1);
        ParameterLabels= Results(1).FitResults(h).ParameterMatrix;
        ResultsMatrix{1,1} = Results(1).FitResults(h).FitType;
        for k = 1:NumberParameters
            if k == NumberParameters -1
                ResultsMatrix{1,k+1} = 'AIC';
            elseif k == NumberParameters
                ResultsMatrix{1,k+1} = 'NLL';
            else
                ResultsMatrix{1,k+1} = ParameterLabels{1,k};
            end
        end

        for b = 1:NumberFiles
            Name = Results(b).Name;
            Name = Name(1:min(MinNameLength,length(Name)));
        %     ResultsMatrix{1,b+1} = Name;
            CurrentParameterMatrix = Results(b).FitResults(h).ParameterMatrix;
            ResultsMatrix{b+1,1} = Name;

            for k = 1:NumberParameters
                if k == NumberParameters -1
                    ResultsMatrix{b+1,k+1} = Results(b).FitResults(h).AIC;
                elseif k == NumberParameters
                    ResultsMatrix{b+1,k+1} = Results(b).FitResults(h).LogLike;
                else
                    ResultsMatrix{b+1,k+1} = CurrentParameterMatrix{2,k};
                end
                
            end
        end
        
        FitParamMatrix(h).ResultsMatrix = ResultsMatrix;
    end

end

function [ResultsMatrix] = Generate_Results_Matrix(Results,UsefulInfo,ResultsOption,MinNameLength)

NumberFiles = length(Results);

if strcmp(ResultsOption,'Fit Stats')
    NumberFitStats = 3;
    NumberFits = UsefulInfo.NumberFitsToPerform;
    ResultsMatrix = cell(NumberFits+1,NumberFitStats+2,NumberFiles);
    
    for b = 1:NumberFiles
        Name = Results(b).Name;
        Name = Name(1:min(MinNameLength,length(Name)));

        ResultsMatrix{1,1,b} = Name;
        ResultsMatrix{1,2,b} = 'Dif_AIC';
        ResultsMatrix{1,3,b} = 'AIC';
        ResultsMatrix{1,4,b} = 'LogLike';
        ResultsMatrix{1,5,b} = 'R2';
        
        for h = 1:NumberFits
            CurrentFitStatsMatrix = Results(b).FitResults(h).FitStatsMatrix;
            FitType = Results(b).FitResults(h).FitType;
            ResultsMatrix{h+1,1,b} = FitType;
            AICList(h) = CurrentFitStatsMatrix{1,1};
            ResultsMatrix(h+1,3:end,b) = CurrentFitStatsMatrix(1,:);
        end
        
        for h = 1:NumberFits
            ResultsMatrix{h+1,2,b} = AICList(h) - min(AICList);
        end
    end
   
elseif strcmp(ResultsOption,'Multiple Fit, AIC')
    
    NumberResultsColumns = UsefulInfo.NumberFitsToPerform;
    ResultsMatrix = cell(NumberFiles +1,NumberResultsColumns + 1);
    ResultsMatrix{1,1} = 'AIC DIFF';
    
    for j = 1:NumberResultsColumns
         ResultsMatrix{1,j+1} = strcat(Results(1).FitResults(j).FitType);
%        ResultsMatrix{1,j+1} = strcat('N = ',num2str(j));
    end
    
    for b = 1:NumberFiles
        Name = Results(b).Name;
        Name = Name(1:min(20,length(Name)));
    %     ResultsMatrix{1,b+1} = Name;
        ResultsMatrix{b+1,1} = Name;

%         for k = 1:NumberResultsColumns
%             ResultsMatrix{b+1,k+1} = Results(b).FitResults(k).AIC;
%         end
            
        for k = 1:NumberResultsColumns
            AIC(k) = Results(b).FitResults(k).AIC;
        end
        for k = 1:NumberResultsColumns
%             ResultsMatrix{b+1,k+1} = AIC(k);
            ResultsMatrix{b+1,k+1} = AIC(k)-min(AIC);
        end
        
    end
    
elseif strcmp(ResultsOption,'Randomness')
    NumberResultsColumns = 2;
    ResultsMatrix = cell(NumberFiles +1,NumberResultsColumns + 1);
    ResultsMatrix{1,1} = '------';
    ResultsMatrix{1,2} = 'RandomParam';
    ResultsMatrix{1,3} = 'Nmin';

    for b = 1:NumberFiles
        Name = Results(b).Name;
        Name = Name(1:min(20,length(Name)));
    %     ResultsMatrix{1,b+1} = Name;
        ResultsMatrix{b+1,1} = Name;

            ResultsMatrix{b+1,2} = Results(b).RandomnessParameter;
            ResultsMatrix{b+1,3} = Results(b).Nmin;
    end
    
elseif strcmp(ResultsOption,'Constrained Best Fit')
    
    NumberResultsColumns = 4;
    NumberofFits = UsefulInfo.NumberFitsToPerform;
    ResultsMatrix = cell(NumberFiles +1,NumberResultsColumns + 1);
    ResultsMatrix{1,1} = 'Constrained Fit';
    ResultsMatrix{1,2} = 'NumberSteps';
    ResultsMatrix{1,3} = 'Tau (s)';
    ResultsMatrix{1,4} = 'Lag Time (s)';
    ResultsMatrix{1,5} = 'R2';

    for b = 1:NumberFiles
        Name = Results(b).Name;
        Name = Name(1:min(20,length(Name)));
    %     ResultsMatrix{1,b+1} = Name;
        ResultsMatrix{b+1,1} = Name;

        AIC = [];
        for k = 1:NumberofFits
            AIC(k) = Results(b).FitResults(k).AIC;
        end
        IndexBestFit = find(AIC == min(AIC));
        
            ResultsMatrix{b+1,2} = Results(b).FitResults(IndexBestFit).NumberSteps;
            ResultsMatrix{b+1,3} = Results(b).FitResults(IndexBestFit).Tau;
            ResultsMatrix{b+1,4} = Results(b).FitResults(IndexBestFit).LagTime;
            ResultsMatrix{b+1,5} = Results(b).FitResults(IndexBestFit).R2;
    end
    
    
end
    


end