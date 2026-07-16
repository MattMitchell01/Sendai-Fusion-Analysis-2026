function [DataToSave, EffectiveDeadTimes] = Combine_Trace_Analysis_Data(DefaultPathname, DataFilenames, NumberDataFiles)
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%

        EffectiveDeadTimes = nan(1,NumberDataFiles);

        for i = 1:NumberDataFiles
            CurrDataFileName = DataFilenames{1,i};
            CurrDataFilePath = strcat(DefaultPathname,CurrDataFileName);

            InputData = open(CurrDataFilePath);

            if isfield(InputData.DataToSave,'OtherDataToSave') && ...
                    isfield(InputData.DataToSave.OtherDataToSave,'EffectiveDeadTime')
                EffectiveDeadTimes(i) = InputData.DataToSave.OtherDataToSave.EffectiveDeadTime;
            end

            if i == 1
                DataToSave = InputData.DataToSave;
            else
                NewDataToAdd = InputData.DataToSave;
                PreviousNumberTraces = length(DataToSave.CombinedAnalyzedTraceData);
                NumbertoAdd = length(NewDataToAdd.CombinedAnalyzedTraceData);


                DataToSave.CombinedAnalyzedTraceData(PreviousNumberTraces+1:NumbertoAdd+PreviousNumberTraces)...
                    = NewDataToAdd.CombinedAnalyzedTraceData;
            end

        end
end

