function [StatsMatrix] = Run_Stats_CDF(AllResults)

NumberFiles = length(AllResults);
StatsMatrix = cell(NumberFiles +1);
StatsMatrix{1,1} = 'KS2 Matrix';

for b = 1:NumberFiles
    Name = AllResults(b).CDFData.Name;
    Name = Name(1:min(12,length(Name)));
    StatsMatrix{1,b+1} = Name;
    StatsMatrix{b+1,1} = Name;
    Data1 = AllResults(b).CDFData.SortedpHtoFList;
    for k = 1:NumberFiles
        Data2 = AllResults(k).CDFData.SortedpHtoFList;
        AlphaValue = 0.05;
        
        % Check to make sure each set of data is a vector (this is to deal 
        % with cases where we are looking at the total video intensity instead of a CDF)
        if isvector(Data1) && isvector(Data2)
            [TestPassed,PValue] = kstest2(Data1,Data2,'Alpha',AlphaValue);

            if PValue < 1e-6
                PValue = 0;
            end
            StatsMatrix{b+1,k+1} = PValue;        
        else
            StatsMatrix{b+1,k+1} = NaN;
        end
        
    end    

%     FileNumber = b;
%     if FileNumber ~= 1
%         CurrentData = CDFData(FileNumber).SortedpHtoFList;
%         PreviousData = CDFData(FileNumber - 1).SortedpHtoFList;
% 
%         AlphaValue = 0.05;
%         [TestPassed,PValue] = kstest2(CurrentData,PreviousData,'Alpha',AlphaValue);
% 
%         if TestPassed
%             Results(FileNumber).KS2TestWPrev =  strcat('Diff, p < ',num2str(AlphaValue));
%         else 
%             Results(FileNumber).KS2TestWPrev =  strcat('Same, p > ',num2str(AlphaValue));
%         end
% 
%         Results(FileNumber).KS2PValue = PValue;
%     else 
%         Results(FileNumber).KS2TestWPrev = 'N/A';
%         Results(FileNumber).KS2PValue = NaN;
%     end

end


end