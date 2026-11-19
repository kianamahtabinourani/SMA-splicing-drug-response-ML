%% 06 Define Biological Classes
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% SMN2: Define biological classes (drug-responsive / non-responsive, weak / strong splice site)
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% Input / Output paths
inFile = fullfile(projectRoot, "Merged data", "SMN2_full_splice_regulatory_data.xlsx");

outDir = fullfile(projectRoot, "Biological class");
outFile = fullfile(outDir, "SMN2_biological_class.xlsx");

if ~exist(outDir, "dir")
    mkdir(outDir);
end

%% Read data
T = readtable(inFile);

%% ---- Ensure needed columns exist (adjust names here if your columns differ) ----
% Drug response columns (from your PSI+stats pipeline)
dpsiRG  = T.dPSI_RG_mean;
dpsiNVS = T.dPSI_NVS_mean;

fdrRG  = T.FDR_RG_vs_DMSO;
fdrNVS = T.FDR_NVS_vs_DMSO;

% Baseline PSI under DMSO
psiDMSO = T.PSI_DMSO_mean;

% Splice-site sequence features (from your splice-regulatory table)
hasGT   = T.has_GT;                 % logical (true/false)
mmCount = T.donor_mismatch_count;   % integer
gcDonor = T.donor10_GC;             % GC% of donor10

%% =========================
% 1) Drug responsiveness labels
% =========================
% Thresholds (simple + standard; you can change them later easily)
DELTA_STRONG = 0.10;   % biologically meaningful change
DELTA_WEAK   = 0.05;   % weaker but noticeable change

FDR_STRICT = 0.05;
FDR_LOOSE  = 0.10;

% RG label
RG_label = strings(height(T),1);
for i = 1:height(T)
    if abs(dpsiRG(i)) >= DELTA_STRONG && fdrRG(i) < FDR_STRICT
        RG_label(i) = "RG_responsive";
    elseif abs(dpsiRG(i)) >= DELTA_WEAK && fdrRG(i) < FDR_LOOSE
        RG_label(i) = "RG_weakly_responsive";
    else
        RG_label(i) = "RG_non_responsive";
    end
end

% NVS label
NVS_label = strings(height(T),1);
for i = 1:height(T)
    if abs(dpsiNVS(i)) >= DELTA_STRONG && fdrNVS(i) < FDR_STRICT
        NVS_label(i) = "NVS_responsive";
    elseif abs(dpsiNVS(i)) >= DELTA_WEAK && fdrNVS(i) < FDR_LOOSE
        NVS_label(i) = "NVS_weakly_responsive";
    else
        NVS_label(i) = "NVS_non_responsive";
    end
end

%% =========================
% 2) Splice site strength labels (weak / strong / intermediate)
% =========================
% Simple biologically interpretable rules using:
% - baseline PSI (DMSO)
% - canonical GT
% - mismatch count
% - donor GC (proxy for local composition)

spliceStrength = strings(height(T),1);

for i = 1:height(T)

    % Strong splice site criteria
    isStrong = (psiDMSO(i) >= 0.70) && (hasGT(i) == true) && (mmCount(i) <= 1);

    % Weak splice site criteria
    isWeak = (psiDMSO(i) <= 0.30) || (hasGT(i) == false) || (mmCount(i) >= 3);

    % Optional: use donor GC to support weak classification
    % (extremely high or low GC can indicate unusual donor context)
    if gcDonor(i) < 20 || gcDonor(i) > 80
        isWeak = true;
    end

    if isStrong
        spliceStrength(i) = "strong_splice_site";
    elseif isWeak
        spliceStrength(i) = "weak_splice_site";
    else
        spliceStrength(i) = "intermediate_splice_site";
    end
end

%% Add label columns to table
T.RG_drug_response_class  = RG_label;
T.NVS_drug_response_class = NVS_label;
T.splice_site_strength_class = spliceStrength;

%% Save
writetable(T, outFile);
disp("Saved biological class table to:");
disp(outFile);