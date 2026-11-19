%% 07 Filter Branaplam Classes
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% Extract NVS (Branaplam) biological classes – categorical version
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% Paths
inFile = fullfile(projectRoot, "Biological class", "SMN2_biological_class.xlsx");
outDir = fullfile(projectRoot, "Biological class");
outFile = fullfile(outDir, "Branaplam_biological_class.xlsx");

%% Read table
T = readtable(inFile);

%% Extract required columns
ss        = T.ss;
NVS_class = string(T.NVS_drug_response_class);      % categorical labels
strength  = string(T.splice_site_strength_class);

%% Define NVS non-responsive flag (kept for convenience)
NVS_non_responsive = NVS_class == "NVS_non_responsive";

%% Build output table (4 columns)
Out = table();
Out.ss = ss;
Out.NVS_response_class = NVS_class;          % keeps non / weak / responsive
Out.NVS_non_responsive = NVS_non_responsive; % logical flag
Out.splice_site_strength = strength;

%% Save
writetable(Out, outFile);
disp("Saved NVS biological class file to:");
disp(outFile);