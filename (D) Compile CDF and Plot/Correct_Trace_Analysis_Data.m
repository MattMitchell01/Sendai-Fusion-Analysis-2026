function [DataToSave] = Correct_Trace_Analysis_Data(DataToSave, Options)
%
% ------------------------------------------------------------------------
% Modified by Matthew D. Mitchell, Rawle Lab, Williams College, 2026
% (parameter / configuration updates to Bob Rawle's original script;
% original pipeline: Rawle et al., Biophysical Journal 2016,
% doi:10.1016/j.bpj.2016.05.048).
% ------------------------------------------------------------------------
%

    NumbertoCheck = length(DataToSave.CombinedAnalyzedTraceData);
    for k = 1:NumbertoCheck
        CurrentFusionData = DataToSave.CombinedAnalyzedTraceData(k).FusionData;

        if strcmp(CurrentFusionData.Designation,'1 Fuse') && ...
                CurrentFusionData.BindtoFusionTime(1) < Options.HiEndOfExclusion && ...
                CurrentFusionData.BindtoFusionTime(1) > Options.LowEndOfExclusion 

            DataToSave.CombinedAnalyzedTraceData(k).Designation = 'User Erased';
            DataToSave.CombinedAnalyzedTraceData(k).FusionData.Designation = ...
                'User Erased';
        end
    end

end
