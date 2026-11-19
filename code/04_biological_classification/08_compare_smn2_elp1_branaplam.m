%% 08 Compare Smn2 Elp1 Branaplam
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% Build NVS (Branaplam) contrast groups: SMN2 vs ELP1
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% Input files (confirmed correct)
smn2File = fullfile(projectRoot, "Biological class", "SMN2_biological_class.xlsx");
elp1File = fullfile(projectRoot, "Biological class", "ELP1_biological_class.xlsx");

%% Output folders (same as before)
outDir1 = fullfile(projectRoot, "Drug effect");
outDir2 = fullfile(projectRoot, "off target susceptibility");

if ~exist(outDir1, "dir"), mkdir(outDir1); end
if ~exist(outDir2, "dir"), mkdir(outDir2); end

%% Read tables
SMN2 = readtable(smn2File);
ELP1 = readtable(elp1File);

%% Ensure key column exists
if ~any(strcmp(SMN2.Properties.VariableNames, 'ss'))
    SMN2.Properties.VariableNames{1} = 'ss';
end
if ~any(strcmp(ELP1.Properties.VariableNames, 'ss'))
    ELP1.Properties.VariableNames{1} = 'ss';
end

%% Convert relevant columns to string (safe comparisons)
SMN2.NVS_drug_response_class     = string(SMN2.NVS_drug_response_class);
SMN2.splice_site_strength_class  = string(SMN2.splice_site_strength_class);

ELP1.NVS_drug_response_class     = string(ELP1.NVS_drug_response_class);
ELP1.splice_site_strength_class  = string(ELP1.splice_site_strength_class);

%% ============================================================
% GROUP DEFINITIONS (Branaplam / NVS contrasts)
% ============================================================

% --- SMN2 group: strong splice site strength AND NVS non-responsive
SMN2_group = SMN2( ...
    SMN2.splice_site_strength_class == "strong_splice_site" & ...
    SMN2.NVS_drug_response_class    == "NVS_non_responsive", :);

% --- ELP1 group: weak splice site strength AND NVS responsive or weakly responsive
% Exclude NVS_not_testable explicitly
ELP1_group = ELP1( ...
    ELP1.splice_site_strength_class == "weak_splice_site" & ...
    (ELP1.NVS_drug_response_class == "NVS_responsive" | ...
     ELP1.NVS_drug_response_class == "NVS_weakly_responsive") & ...
    ELP1.NVS_drug_response_class ~= "NVS_not_testable", :);

%% Keep only columns needed for downstream sequence analysis
keepCols = ["ss","NVS_drug_response_class","splice_site_strength_class"];

SMN2_group = SMN2_group(:, intersect(keepCols, string(SMN2_group.Properties.VariableNames), 'stable'));
ELP1_group = ELP1_group(:, intersect(keepCols, string(ELP1_group.Properties.VariableNames), 'stable'));

%% Save outputs (same groups saved into both requested directories)
smn2OutName = "SMN2_NVS_strong_nonresponsive.xlsx";
elp1OutName = "ELP1_NVS_weak_responsive.xlsx";

% Contrast set saved to "Drug effect"
writetable(SMN2_group, fullfile(outDir1, smn2OutName));
writetable(ELP1_group, fullfile(outDir1, elp1OutName));

% Contrast set saved to "off target susceptibility"
writetable(SMN2_group, fullfile(outDir2, smn2OutName));
writetable(ELP1_group, fullfile(outDir2, elp1OutName));

%% Report counts (sanity check)
fprintf("Saved NVS contrast groups.\n");
fprintf("SMN2 strong & NVS non-responsive: %d rows\n", height(SMN2_group));
fprintf("ELP1 weak & NVS responsive/weakly responsive: %d rows\n", height(ELP1_group));

disp("Outputs saved to:");
disp(outDir1);
disp(outDir2);