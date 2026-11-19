%% 02 Extract Psi Features
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% SMN2: Extract PSI-derived features (including replicate-level deltaPSI)
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% Paths
inCsv  = fullfile(projectRoot, "data", "raw", "GSE221868_smn2_psi.csv");
outDir = fullfile(projectRoot, "PSI features");
outXlsx = fullfile(outDir, "SMN2_PSI_features.xlsx");

% Create output folder if it doesn't exist
if ~exist(outDir, "dir")
    mkdir(outDir);
end

%% Read table
T = readtable(inCsv);

% First column = 5' splice site sequences
T.Properties.VariableNames{1} = 'ss';
ss = upper(string(T.ss));

vars = string(T.Properties.VariableNames);

%% Select columns for each condition by prefix
dmsoMask = startsWith(vars, "smn2_dmso_");
rgMask   = startsWith(vars, "smn2_rg_");    % risdiplam (RG)
nvsMask  = startsWith(vars, "smn2_nvs_");   % branaplam (NVS)

if sum(dmsoMask)==0 || sum(rgMask)==0 || sum(nvsMask)==0
    error("Could not find expected column prefixes (smn2_dmso_, smn2_rg_, smn2_nvs_).");
end

D = table2array(T(:, dmsoMask));  % 285 x 9
R = table2array(T(:, rgMask));    % 285 x 9
N = table2array(T(:, nvsMask));   % 285 x 9

% Safety: must be paired (same number of columns)
if size(D,2) ~= size(R,2) || size(D,2) ~= size(N,2)
    error("Mismatch in number of replicate columns: DMSO=%d, RG=%d, NVS=%d", ...
        size(D,2), size(R,2), size(N,2));
end

%% Helper for safe division
eps0 = 1e-9;

%% Mean PSI and SD PSI per condition
PSI_DMSO_mean = mean(D, 2, "omitnan");
PSI_RG_mean   = mean(R, 2, "omitnan");
PSI_NVS_mean  = mean(N, 2, "omitnan");

PSI_DMSO_sd = std(D, 0, 2, "omitnan");
PSI_RG_sd   = std(R, 0, 2, "omitnan");
PSI_NVS_sd  = std(N, 0, 2, "omitnan");

%% Delta PSI (paired per (library,replicate))
% Each column is assumed to correspond to matched (lib,rep) ordering.
dRG  = R - D;   % 285 x 9
dNVS = N - D;   % 285 x 9

% Summary delta PSI
dPSI_RG_mean  = mean(dRG,  2, "omitnan");
dPSI_NVS_mean = mean(dNVS, 2, "omitnan");

dPSI_RG_sd  = std(dRG,  0, 2, "omitnan");
dPSI_NVS_sd = std(dNVS, 0, 2, "omitnan");

%% Replicate-level delta PSI features (explicit)
dPSI_RG_reps  = array2table(dRG);
dPSI_NVS_reps = array2table(dNVS);

for i = 1:width(dPSI_RG_reps)
    dPSI_RG_reps.Properties.VariableNames{i}  = sprintf('dPSI_RG_rep%d', i);
    dPSI_NVS_reps.Properties.VariableNames{i} = sprintf('dPSI_NVS_rep%d', i);

end

%% CV per condition (SD/mean)
CV_DMSO = PSI_DMSO_sd ./ max(abs(PSI_DMSO_mean), eps0);
CV_RG   = PSI_RG_sd   ./ max(abs(PSI_RG_mean),   eps0);
CV_NVS  = PSI_NVS_sd  ./ max(abs(PSI_NVS_mean),  eps0);

%% SD/CV per condition
SD_over_CV_DMSO = PSI_DMSO_sd ./ max(CV_DMSO, eps0);
SD_over_CV_RG   = PSI_RG_sd   ./ max(CV_RG,   eps0);
SD_over_CV_NVS  = PSI_NVS_sd  ./ max(CV_NVS,  eps0);

%% Build output table
Out = table();
Out.ss = ss;

% Mean PSI
Out.PSI_DMSO_mean = PSI_DMSO_mean;
Out.PSI_RG_mean   = PSI_RG_mean;
Out.PSI_NVS_mean  = PSI_NVS_mean;

% SD PSI
Out.PSI_DMSO_sd = PSI_DMSO_sd;
Out.PSI_RG_sd   = PSI_RG_sd;
Out.PSI_NVS_sd  = PSI_NVS_sd;

% Delta PSI (summary)
Out.dPSI_RG_mean  = dPSI_RG_mean;
Out.dPSI_NVS_mean = dPSI_NVS_mean;
Out.dPSI_RG_sd    = dPSI_RG_sd;
Out.dPSI_NVS_sd   = dPSI_NVS_sd;

% CV
Out.CV_DMSO = CV_DMSO;
Out.CV_RG   = CV_RG;
Out.CV_NVS  = CV_NVS;

% SD/CV
Out.SD_over_CV_DMSO = SD_over_CV_DMSO;
Out.SD_over_CV_RG   = SD_over_CV_RG;
Out.SD_over_CV_NVS  = SD_over_CV_NVS;

% Add replicate-level delta PSI columns
Out = [Out, dPSI_RG_reps, dPSI_NVS_reps];

%% Save
writetable(Out, outXlsx);
disp("Saved: " + outXlsx);