function [InputData,CurrentFilename] = Get_Data_From_Other_Source(SourceNumber)

if SourceNumber == 1
    DataFromHistogram = [0.01, 0.048, 0.05, 0.0325, 0.018, 0.009, 0.012, 0.004, 0.002, 0.005, 0.002, 0.001, 0.001, 0.0025, 0.0001, 0.0001, 0.0025];
    TimePoints = 5:5:length(DataFromHistogram)*5;
    NumberVirus = 309*1/sum(DataFromHistogram);
    DataSource = 'Floyd, 2008, Fig 3A, pH 4,6.mat';
elseif SourceNumber == 2
    DataFromHistogram = 0.001*[0.4, 1.8, 4.2, 5.5, 8.2, 8.3, 6.5, 5.9, 5.9, 5.2, 6.8, 5.7, 2.5, 4.0, 1.6, 0.2, 0.5, 1.6, 1.0];
    TimeSpacing = 12.5;
    TimePoints = TimeSpacing:TimeSpacing:length(DataFromHistogram)*TimeSpacing;
    NumberVirus = 380;
    NumberVirus = NumberVirus*1/sum(DataFromHistogram);
    DataSource = 'Ivanovic,2013,Fig1C, pH 5,5.mat';
end

CurrentFilename = DataSource;

% compile into data format needed for fitting program
NumberRecorded = 0;
for y= 1:length(DataFromHistogram)
    CurrentNumberPoints = round(NumberVirus*DataFromHistogram(y));
    CurrentTimePoint = TimePoints(y);
    
    if CurrentNumberPoints >= 1
        for n= 1:CurrentNumberPoints
            NumberRecorded = NumberRecorded + 1;
            InputData.Useful_Data_To_Save.Combined_Fuse1_Data(NumberRecorded).pHtoFuse1Time = CurrentTimePoint;
        end
    end
end

end