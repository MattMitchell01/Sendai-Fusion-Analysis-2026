% This is a quick code to extract some CDF data into a separate file for a specific purpose

for w = 1:NumberDataFiles
    CDFData(w) = AllResults(w).CDFData;
end

for w = 1:NumberDataFiles
    CDFData(w).pH = [];
    CDFData(w).EfficiencyAfter = ResultsReport(w).PercentFuse1*100;
    CDFData(w).EfficiencyBefore = [];
end

save('ZikaDataToFit.mat','CDFData')

TestData = CDFData(1).SortedpHtoFList;
Points = min(TestData):.001:max(TestData);
% Test2 = NaN(5000, 1);
% TestData = [TestData; Test2];
[f,xi] = ksdensity(TestData,Points,'Support','positive');
figure
plot(xi,f);
trapz(xi,f)