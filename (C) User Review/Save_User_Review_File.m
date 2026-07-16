function Save_User_Review_File(UserReviewFilePath, DataToSave, ReviewQueue, CurrentQueuePosition, ReviewMeta)
%
% ------------------------------------------------------------------------
% Written by Matthew D. Mitchell, Rawle Lab, Williams College, 2026.
% New script developed as part of the 2026 overhaul and expansion of the
% Sendai fusion analysis pipeline originally created by Bob Rawle
% (Kasson Lab, University of Virginia, 2016; Rawle et al., Disentangling
% Viral Membrane Fusion from Receptor Binding Using Synthetic DNA-Lipid
% Conjugates, Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% pipeline subsequently updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% Save_User_Review_File  Single place that writes the User Review File.
%
% Used both for the first-run bootstrap and for every later round's
% save-in-place, so the variable contract can never drift between them.

    UserReviewDir = fileparts(UserReviewFilePath);
    if ~isempty(UserReviewDir) && exist(UserReviewDir, 'dir') == 0
        mkdir(UserReviewDir);
    end

    save(UserReviewFilePath, 'DataToSave', 'ReviewQueue', 'CurrentQueuePosition', 'ReviewMeta');
end
