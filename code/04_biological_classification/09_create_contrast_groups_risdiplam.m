%% 09 Create Contrast Groups Risdiplam
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% RG (Risdiplam): Build contrast groups for SMN2 vs SMN2 and SMN2 vs ELP1
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% Input files
smn2File = fullfile(projectRoot, "Biological class", "SMN2_biological_class.xlsx");
elp1File = fullfile(projectRoot, "Biological class", "ELP1_biological_class.xlsx");  % (ELP1 file; name doesn't matter)

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

%% Make sure RG class exists
if ~any(strcmp(SMN2.Properties.VariableNames, "RG_drug_response_class"))
    error("SMN2 file does not contain RG_drug_response_class.");
end
if ~any(strcmp(ELP1.Properties.VariableNames, "RG_drug_response_class"))
    error("ELP1 file does not contain RG_drug_response_class.");
end

%% Convert to string for safe comparisons
SMN2.RG_drug_response_class = string(SMN2.RG_drug_response_class);
ELP1.RG_drug_response_class = string(ELP1.RG_drug_response_class);

%% Define responsive / non-responsive masks
SMN2_resp = (SMN2.RG_drug_response_class == "RG_responsive") | ...
            (SMN2.RG_drug_response_class == "RG_weakly_responsive");

SMN2_non  = (SMN2.RG_drug_response_class == "RG_non_responsive");

ELP1_non  = (ELP1.RG_drug_response_class == "RG_non_responsive");  % (if you also have RG_not_testable, it won't be included)

%% Build groups
SMN2_responsive_group    = SMN2(SMN2_resp, :);
SMN2_nonresponsive_group = SMN2(SMN2_non,  :);

ELP1_nonresponsive_group = ELP1(ELP1_non,  :);

%% Save outputs
writetable(SMN2_responsive_group,    fullfile(outDir, "RG_SMN2_responsive.xlsx"));
writetable(SMN2_nonresponsive_group, fullfile(outDir, "RG_SMN2_nonresponsive.xlsx"));
writetable(ELP1_nonresponsive_group, fullfile(outDir, "RG_ELP1_nonresponsive.xlsx"));

% Also save the cross-gene contrast as a paired “set” (two tabs would require writecell/write to sheets;
% simplest is two files, but we print counts for clarity.)
fprintf("Saved RG contrast groups to: %s\n", outDir);
fprintf("SMN2 responsive: %d\n", height(SMN2_responsive_group));
fprintf("SMN2 non-responsive: %d\n", height(SMN2_nonresponsive_group));
fprintf("ELP1 non-responsive: %d\n", height(ELP1_nonresponsive_group));
