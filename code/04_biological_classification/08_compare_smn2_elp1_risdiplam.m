%% 08 Compare Smn2 Elp1 Risdiplam
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% Build RG (Risdiplam) contrast groups: SMN2 vs ELP1
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% Input files (confirmed correct)
smn2File = fullfile(projectRoot, "Biological class", "SMN2_biological_class.xlsx");
elp1File = fullfile(projectRoot, "Biological class", "ELP1_biological_class.xlsx");

%% Output folders (as requested)
outDir1 = fullfile(projectRoot, "Drug effect");
outDir2 = fullfile(projectRoot, "off target susceptibility");

if ~exist(outDir1, "dir"), mkdir(outDir1); end
if ~exist(outDir2, "dir"), mkdir(outDir2); end

%% Read tables
SMN2 = readtable(smn2File);
ELP1 = readtable(elp1File);

%% Ensure key columns exist (rename first column to 'ss' if needed)
if ~any(strcmp(SMN2.Properties.VariableNames, 'ss'))
    SMN2.Properties.VariableNames{1} = 'ss';
end
if ~any(strcmp(ELP1.Properties.VariableNames, 'ss'))
    ELP1.Properties.VariableNames{1} = 'ss';
end

%% Pull relevant columns as strings (safer for comparisons)
SMN2.RG_drug_response_class       = string(SMN2.RG_drug_response_class);
SMN2.splice_site_strength_class   = string(SMN2.splice_site_strength_class);

ELP1.RG_drug_response_class       = string(ELP1.RG_drug_response_class);
ELP1.splice_site_strength_class   = string(ELP1.splice_site_strength_class);

%% ============================================================
% GROUP DEFINITIONS (your requested biological contrasts)
% ============================================================

% --- SMN2 group: strong splice site strength AND RG non-responsive
SMN2_group = SMN2( ...
    SMN2.splice_site_strength_class == "strong_splice_site" & ...
    SMN2.RG_drug_response_class     == "RG_non_responsive", :);

% --- ELP1 group: weak splice site strength AND RG responsive
% Exclude RG_not_testable explicitly
ELP1_group = ELP1( ...
    ELP1.splice_site_strength_class == "weak_splice_site" & ...
    (ELP1.RG_drug_response_class == "RG_responsive" | ...
     ELP1.RG_drug_response_class == "RG_weakly_responsive") & ...
    ELP1.RG_drug_response_class ~= "RG_not_testable", :);

%% Keep only the columns needed for downstream sequence analysis
keepCols = ["ss","RG_drug_response_class","splice_site_strength_class"];

SMN2_group = SMN2_group(:, intersect(keepCols, string(SMN2_group.Properties.VariableNames), 'stable'));
ELP1_group = ELP1_group(:, intersect(keepCols, string(ELP1_group.Properties.VariableNames), 'stable'));

%% Save outputs (same groups saved into both requested directories)
smn2OutName = "SMN2_RG_strong_nonresponsive.xlsx";
elp1OutName = "ELP1_RG_weak_responsive.xlsx";

% Contrast set saved to "Drug effect"
writetable(SMN2_group, fullfile(outDir1, smn2OutName));
writetable(ELP1_group, fullfile(outDir1, elp1OutName));

% Contrast set saved to "off target susceptibility"
writetable(SMN2_group, fullfile(outDir2, smn2OutName));
writetable(ELP1_group, fullfile(outDir2, elp1OutName));

%% Report counts (quick sanity check)
fprintf("Saved RG contrast groups.\n");
fprintf("SMN2 strong & RG non-responsive: %d rows\n", height(SMN2_group));
fprintf("ELP1 weak & RG responsive (not_testable excluded): %d rows\n", height(ELP1_group));

disp("Outputs saved to:");
disp(outDir1);
disp(outDir2);