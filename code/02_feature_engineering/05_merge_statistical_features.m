%% 05 Merge Statistical Features
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% Merge SMN2 merged features with statistical features
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% Paths
mergedFile = fullfile(projectRoot, "Merged data", "SMN2_merged.xlsx");
statFile   = fullfile(projectRoot, "Statistical features", "SMN2_statistical_features.xlsx");

outDir  = fullfile(projectRoot, "Merged data");
outFile = fullfile(outDir, "SMN2_full_splice_regulatory_data.xlsx");

%% Read tables
T_main = readtable(mergedFile);
T_stat = readtable(statFile);

%% Ensure 'ss' is the key column in both tables
if ~any(strcmp(T_main.Properties.VariableNames, 'ss'))
    T_main.Properties.VariableNames{1} = 'ss';
end

if ~any(strcmp(T_stat.Properties.VariableNames, 'ss'))
    T_stat.Properties.VariableNames{1} = 'ss';
end

%% Merge (inner join keeps only matching splice sites)
T_full = innerjoin(T_main, T_stat, 'Keys', 'ss');

%% Save
writetable(T_full, outFile);
disp("Final merged dataset saved to:");
disp(outFile);