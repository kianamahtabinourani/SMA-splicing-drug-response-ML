%% 03 Merge Psi And Regulatory Features
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% Merge SMN2 splice-regulatory features with PSI-derived features
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% Paths
spliceFile = fullfile(projectRoot, "Full splice regulatory dataset", "SMN2_splice_regulatory_features_ONLY.xlsx");
psiFile    = fullfile(projectRoot, "PSI features", "SMN2_PSI_features.xlsx");

outDir  = fullfile(projectRoot, "Merged data");
outFile = fullfile(outDir, "SMN2_merged.xlsx");

% Create output folder if needed
if ~exist(outDir, "dir")
    mkdir(outDir);
end

%% Read tables
T_splice = readtable(spliceFile);
T_psi    = readtable(psiFile);

%% Ensure 'ss' is the key in both tables
% (Rename first column if needed)
if ~any(strcmp(T_splice.Properties.VariableNames, 'ss'))
    T_splice.Properties.VariableNames{1} = 'ss';
end
if ~any(strcmp(T_psi.Properties.VariableNames, 'ss'))
    T_psi.Properties.VariableNames{1} = 'ss';
end

%% Merge (inner join keeps only matching ss values)
T_merged = innerjoin(T_splice, T_psi, 'Keys', 'ss');

%% Save
writetable(T_merged, outFile);
disp("Merged file saved to: " + outFile);