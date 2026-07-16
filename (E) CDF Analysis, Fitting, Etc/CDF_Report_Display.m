function CDF_Report_Display(CDFReports)

    NumberDataFiles = length(CDFReports);

    disp('   --------------CDF Reports---------------')
    for FileNumber = 1:NumberDataFiles
        disp(CDFReports(FileNumber))
    end

    disp('   ----------------------------------------')
end