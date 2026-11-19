%% 04 Compute Statistical Features
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% SMN2: Wilcoxon signed-rank p-values + FDR (RG vs DMSO, NVS vs DMSO)
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% Paths
inCsv  = fullfile(projectRoot, "data", "raw", "GSE221868_smn2_psi.csv");
outDir = fullfile(projectRoot, "Statistical features");
outXlsx = fullfile(outDir, "SMN2_statistical_features.xlsx");

if ~exist(outDir, "dir")
    mkdir(outDir);
end

%% Read SMN2 PSI table
T = readtable(inCsv);
T.Properties.VariableNames{1} = 'ss';   % rename first column safely
ss = upper(string(T.ss));

vars = string(T.Properties.VariableNames);

%% Identify columns by prefix (matches your file naming)
dmsoMask = startsWith(vars, "smn2_dmso_");
rgMask   = startsWith(vars, "smn2_rg_");
nvsMask  = startsWith(vars, "smn2_nvs_");

if sum(dmsoMask)==0 || sum(rgMask)==0 || sum(nvsMask)==0
    error("Could not find expected columns: smn2_dmso_, smn2_rg_, smn2_nvs_.");
end

D = table2array(T(:, dmsoMask));  % 285 x 9
R = table2array(T(:, rgMask));    % 285 x 9
N = table2array(T(:, nvsMask));   % 285 x 9

% Ensure paired columns count matches (required for signed-rank pairing)
if size(D,2) ~= size(R,2) || size(D,2) ~= size(N,2)
    error("Mismatch in paired columns: DMSO=%d, RG=%d, NVS=%d", size(D,2), size(R,2), size(N,2));
end

%% Wilcoxon signed-rank tests (two-sided), per splice site
n = size(D,1);
p_RG  = nan(n,1);
p_NVS = nan(n,1);

for i = 1:n
    dRG  = R(i,:) - D(i,:);   % paired differences (RG - DMSO)
    dNVS = N(i,:) - D(i,:);   % paired differences (NVS - DMSO)

    % Remove NaNs if any (keep pairing by working on differences)
    dRG  = dRG(~isnan(dRG));
    dNVS = dNVS(~isnan(dNVS));

    % If all diffs are zero or too few points, signrank may return NaN; handle gracefully
    if numel(dRG) >= 3 && any(dRG ~= 0)
        p_RG(i) = signrank(dRG, 0, 'tail', 'both');   % test median(diff)=0
    else
        p_RG(i) = NaN;
    end

    if numel(dNVS) >= 3 && any(dNVS ~= 0)
        p_NVS(i) = signrank(dNVS, 0, 'tail', 'both');
    else
        p_NVS(i) = NaN;
    end
end

%% FDR (Benjamini–Hochberg) for each drug comparison
q_RG  = bh_fdr(p_RG);
q_NVS = bh_fdr(p_NVS);

%% Build output table
Out = table();
Out.ss = ss;

Out.p_RG_vs_DMSO  = p_RG;
Out.FDR_RG_vs_DMSO = q_RG;

Out.p_NVS_vs_DMSO  = p_NVS;
Out.FDR_NVS_vs_DMSO = q_NVS;

%% Save
writetable(Out, outXlsx);
disp("Saved: " + outXlsx);

%% -------- helper: BH-FDR --------
function q = bh_fdr(p)
%BH_FDR Benjamini-Hochberg FDR (q-values) for a vector of p-values.
% Handles NaNs by leaving q as NaN in those positions.

    p = p(:);
    q = nan(size(p));

    ok = ~isnan(p);
    p0 = p(ok);

    if isempty(p0)
        return;
    end

    [ps, idx] = sort(p0, 'ascend');
    m = numel(ps);

    % BH: q_i = p_i * m / i, then enforce monotonicity from the end
    q_sorted = ps .* (m ./ (1:m)');
    q_sorted = min(q_sorted, 1);

    for i = m-1:-1:1
        q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
    end

    % Put back in original order
    q0 = nan(size(p0));
    q0(idx) = q_sorted;

    q(ok) = q0;
end