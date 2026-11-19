%% 07 Filter Risdiplam Classes
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% Extract RG (Risdiplam) biological classes – categorical version
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% Paths
inFile = fullfile(projectRoot, "Biological class", "SMN2_biological_class.xlsx");
outDir = fullfile(projectRoot, "Biological class");
outFile = fullfile(outDir, "Risdipla_biological_class.xlsx");

%% Read table
T = readtable(inFile);

%% Extract required columns
ss        = T.ss;
RG_class  = string(T.RG_drug_response_class);   % categorical labels
strength  = string(T.splice_site_strength_class);

%% Define RG non-responsive flag (kept for convenience)
RG_non_responsive = RG_class == "RG_non_responsive";

%% Build output table (4 columns)
Out = table();
Out.ss = ss;
Out.RG_response_class = RG_class;      % keeps non / weak / responsive
Out.RG_non_responsive = RG_non_responsive;
Out.splice_site_strength = strength;

%% Save
writetable(Out, outFile);
disp("Saved RG biological class file to:");
disp(outFile);