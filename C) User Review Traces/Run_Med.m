function [CurrTrace_Corrected] = Run_Med(CurrTrace_Corrected,Options)
    RunMedHalfLength = Options.RunMedHalfLength;
    TraceLength = length(CurrTrace_Corrected);
    OldTrace = CurrTrace_Corrected;
    NewTrace = zeros(size(CurrTrace_Corrected));

    for n = 1:TraceLength
        if n - RunMedHalfLength < 1
            NewTrace(n) = median(OldTrace(1:n+RunMedHalfLength));
        elseif n + RunMedHalfLength > TraceLength
            NewTrace(n) = median(OldTrace(n-RunMedHalfLength:end));
        else
            NewTrace(n) = median(OldTrace(n-RunMedHalfLength:n+RunMedHalfLength));
        end
        
    end

    CurrTrace_Corrected = NewTrace;

end




