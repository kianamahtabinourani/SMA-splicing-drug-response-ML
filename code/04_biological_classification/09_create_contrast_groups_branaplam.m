%% 09 Create Contrast Groups Branaplam
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% NVS (Branaplam): Build contrast groups for SMN2 vs SMN2 and SMN2 vs ELP1
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% Input files (same directories)
smn2File = fullfile(projectRoot, "Biological class", "SMN2_biological_class.xlsx");
elp1File = fullfile(projectRoot, "Biological class", "ELP1_biological_class.xlsx");

%% Output folder
outDir = fullfile(projectRoot, "Contrast group");
if ~exist(outDir, "dir"), mkdir(outDir); end

%% Read tables
SMN2 = readtable(smn2File);
ELP1 = readtable(elp1File);

%% Ensure key column 'ss' exists
if ~any(strcmp(SMN2.Properties.VariableNames, 'ss'))
    SMN2.Properties.VariableNames{1} = 'ss';
end
if ~any(strcmp(ELP1.Properties.VariableNames, 'ss'))
    ELP1.Properties.VariableNames{1} = 'ss';
end

%% Make sure NVS class exists
if ~any(strcmp(SMN2.Properties.VariableNames, "NVS_drug_response_class"))
    error("SMN2 file does not contain NVS_drug_response_class.");
end
if ~any(strcmp(ELP1.Properties.VariableNames, "NVS_drug_response_class"))
    error("ELP1 file does not contain NVS_drug_response_class.");
end

%% Convert to string for safe comparisons
SMN2.NVS_drug_response_class = string(SMN2.NVS_drug_response_class);
ELP1.NVS_drug_response_class = string(ELP1.NVS_drug_response_class);

%% Define responsive / non-responsive masks
SMN2_resp = (SMN2.NVS_drug_response_class == "NVS_responsive") | ...
            (SMN2.NVS_drug_response_class == "NVS_weakly_responsive");

SMN2_non  = (SMN2.NVS_drug_response_class == "NVS_non_responsive");

% Exclude not_testable explicitly if present (safe even if none exist)
ELP1_non  = (ELP1.NVS_drug_response_class == "NVS_non_responsive");

%% Build groups
SMN2_responsive_group    = SMN2(SMN2_resp, :);
SMN2_nonresponsive_group = SMN2(SMN2_non,  :);
ELP1_nonresponsive_group = ELP1(ELP1_non,  :);

%% Save outputs (requested names)
writetable(ELP1_nonresponsive_group, fullfile(outDir, "NVS_ELP1_nonresponsive.xlsx"));
writetable(SMN2_nonresponsive_group, fullfile(outDir, "NVS_SMN2_nonresponsive.xlsx"));
writetable(SMN2_responsive_group,    fullfile(outDir, "NVS_SMN2_responsive.xlsx"));

%% Report counts
fprintf("Saved NVS contrast groups to: %s\n", outDir);
fprintf("SMN2 responsive: %d\n", height(SMN2_responsive_group));
fprintf("SMN2 non-responsive: %d\n", height(SMN2_nonresponsive_group));
fprintf("ELP1 non-responsive: %d\n", height(ELP1_nonresponsive_group));