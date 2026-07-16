function [BootstrapMedianMatrix] = Run_Bootstrap_Median(AllResults,Options)

NumberFiles = length(AllResults);
BootstrapMedianMatrix = cell(NumberFiles +1);
BootstrapMedianMatrix{1,1} = 'BootMedian Matrix';

for b = 1:NumberFiles
    Name = AllResults(b).CDFData.Name;
    Name = Name(1:min(12,length(Name)));
    BootstrapMedianMatrix{1,b+1} = Name;
    BootstrapMedianMatrix{b+1,1} = Name;
    Data1 = AllResults(b).CDFData.SortedpHtoFList;
    
    for k = 1:NumberFiles
        Data2 = AllResults(k).CDFData.SortedpHtoFList;
%         AlphaValue = 0.05;
        [PValue] = Bootstrap_Median_Test(Data1,Data2,Options);

        if PValue < 1e-6
            PValue = 0;
        end
        
        BootstrapMedianMatrix{b+1,k+1} = PValue;
    end    

end


end

function [PValue] = Bootstrap_Median_Test(Data1,Data2,Options)
    
    NumberBootstraps = Options.NumberBootstraps;
    
    NumberDataPoints1 = length(Data1);
    IndexMatrix1 = randi(NumberDataPoints1,[NumberDataPoints1,NumberBootstraps]);
    BootstrapMatrix1 = Data1(IndexMatrix1);
    BootstrapMedianValues1 = median(BootstrapMatrix1,1);
%     BootstrapMatrix  = sort(BootstrapMatrix);

    NumberDataPoints2 = length(Data2);
    IndexMatrix2 = randi(NumberDataPoints2,[NumberDataPoints2,NumberBootstraps]);
    BootstrapMatrix2 = Data2(IndexMatrix2);
    BootstrapMedianValues2 = median(BootstrapMatrix2,1);
    
    MedianDifference = BootstrapMedianValues1 - BootstrapMedianValues2;
    
    CutoffValue = 0;
    NumberAboveCutoff = length(MedianDifference(MedianDifference > CutoffValue));
    NumberBelowCutoff = length(MedianDifference(MedianDifference <= CutoffValue));
    PValue = NumberBelowCutoff/(NumberAboveCutoff+NumberBelowCutoff);
end