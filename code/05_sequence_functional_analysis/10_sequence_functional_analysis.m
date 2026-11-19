%% 10 Sequence Functional Analysis
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% Sequence Functional Analysis (ESE/ESS enrichment, SR motif overlap, 5'SS controls)
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% ---- Paths ----
contrastFiles = {
    fullfile(projectRoot, "Contrast group", "NVS_ELP1_nonresponsive.xlsx")
    fullfile(projectRoot, "Contrast group", "NVS_SMN2_nonresponsive.xlsx")
    fullfile(projectRoot, "Contrast group", "NVS_SMN2_responsive.xlsx")
    fullfile(projectRoot, "Contrast group", "RG_ELP1_nonresponsive.xlsx")
    fullfile(projectRoot, "Contrast group", "RG_SMN2_nonresponsive.xlsx")
    fullfile(projectRoot, "Contrast group", "RG_SMN2_responsive.xlsx")
};

smn2FeatFile = fullfile(projectRoot, "Biological class", "SMN2_biological_class.xlsx");
elp1FeatFile = fullfile(projectRoot, "Biological class", "ELP1_biological_class.xlsx");

outDir   = fullfile(projectRoot, "Sequence functional analysis");
motifDir = fullfile(outDir, "motifs");

if ~exist(outDir,"dir"), mkdir(outDir); end

%% ---- Motif lists ----
% Real motif sets:
%   ESE: RESCUE-ESE hexamers (238)
%   ESS: FAS-hex3 hexamers
motifs = loadMotifs_RES_And_FAS(motifDir);

%% ---- SR-protein motif scoring (ESEfinder matrices + thresholds) ----
SR = makeESEfinderMatrices();  % you already saved this function

%% ---- Read feature tables ----
SMN2 = standardizeSSColumn(readtable(smn2FeatFile));
ELP1 = standardizeSSColumn(readtable(elp1FeatFile));

%% We scan motifs on exonSeq (best for ESE/ESS/SR)
seqColumnToScan = "exonSeq";

%% ---- Process each contrast file ----
for i = 1:numel(contrastFiles)
    cf = contrastFiles{i};

    C = standardizeSSColumn(readtable(cf));
    [drug, gene, cls] = parseContrastFilename(cf);

    if gene == "SMN2"
        F = SMN2;
    else
        F = ELP1;
    end

    % Match by ss key
    keyC = upper(string(C.ss));
    keyF = upper(string(F.ss));
    keep = ismember(keyF, keyC);
    Sub = F(keep, :);

    % Safety check: ensure exonSeq exists
    if ~ismember(seqColumnToScan, string(Sub.Properties.VariableNames))
        error("Column '%s' not found. Available columns include: %s", ...
            seqColumnToScan, strjoin(string(Sub.Properties.VariableNames), ", "));
    end

    % Clean sequence for scanning
    seqs = upper(string(Sub.(seqColumnToScan)));
    seqs = replace(seqs,"U","T");
    seqs = regexprep(seqs,"[^ACGT]",""); % keep only A/C/G/T

    % ----- ESE/ESS motif counts -----
    motifFeat = computeMotifFeatures_ESE_ESS(seqs, motifs);

    % ----- SR-protein binding via ESEfinder PWMs -----
    [sf2_hits,  sf2_max]   = scanPWM(seqs, SR.SF2_ASF.matrix, SR.SF2_ASF.threshold);
    [sc35_hits, sc35_max]  = scanPWM(seqs, SR.SC35.matrix,   SR.SC35.threshold);
    [srp40_hits, srp40_max]= scanPWM(seqs, SR.SRp40.matrix,  SR.SRp40.threshold);
    [srp55_hits, srp55_max]= scanPWM(seqs, SR.SRp55.matrix,  SR.SRp55.threshold);

    SRfeat = table(sf2_hits, sf2_max, sc35_hits, sc35_max, srp40_hits, srp40_max, srp55_hits, srp55_max);

    % ----- Keep ID + PSI/drug-response + 5'SS controls if present -----
    keepCols = intersect( ...
        ["ss","donor10_ref","donor10","has_GT","donor_mismatch_count","donor10_GC", ...
         "splice_site_strength_class","donor10_GC","donor10_GC", ...
         "exon7_GC","upIntron200_GC","downIntron200_GC","upIntron_last30_pyFrac", ...
         "PSI_DMSO_mean","PSI_RG_mean","PSI_NVS_mean","dPSI_RG_mean","dPSI_NVS_mean", ...
         "p_RG_vs_DMSO","FDR_RG_vs_DMSO","p_NVS_vs_DMSO","FDR_NVS_vs_DMSO", ...
         "RG_drug_response_class","NVS_drug_response_class"], ...
        string(Sub.Properties.VariableNames), "stable");

    % IMPORTANT: do NOT append ctrlFeat separately (avoids duplicate column names)
    Out = [Sub(:, keepCols), motifFeat, SRfeat];

    outName = sprintf("%s_%s_%s_sequence_functional.xlsx", drug, gene, cls);
    outPath = fullfile(outDir, outName);
    writetable(Out, outPath);
    fprintf("Saved: %s (n=%d)\n", outPath, height(Out));
end

%% ---- Enrichment summary: responsive vs nonresponsive within each (drug,gene) ----
saveEnrichmentSummaries(outDir);
fprintf("\nDONE. Outputs are in:\n%s\n", outDir);

%% ================= Local functions =================

function T = standardizeSSColumn(T)
    if ~ismember("ss", string(T.Properties.VariableNames))
        T.Properties.VariableNames{1} = "ss";
    end
end

function motifs = loadMotifs_RES_And_FAS(motifDir)
    motifs = struct();

    eseFile = fullfile(motifDir, "RESCUE_ESE_238.txt");
    essFile = fullfile(motifDir, "FAS_hex3.txt");

    if ~exist(eseFile, "file")
        error("ESE file not found: %s", eseFile);
    end
    if ~exist(essFile, "file")
        error("ESS file not found: %s", essFile);
    end

    motifs.ESE = readHexamerFile(eseFile);
    motifs.ESS = readHexamerFile(essFile);

    fprintf("Loaded motifs: ESE=%d hexamers, ESS=%d hexamers\n", numel(motifs.ESE), numel(motifs.ESS));
end

function hex = readHexamerFile(fp)
    raw = splitlines(string(fileread(fp)));
    raw = strtrim(raw);
    raw(raw=="") = [];

    tok = strings(0,1);
    for i = 1:numel(raw)
        line = raw(i);

        % skip header/comment-like lines
        if startsWith(line, "#") || startsWith(line, ">") || startsWith(lower(line), "motif")
            continue;
        end

        parts = split(line);
        first = upper(parts(1));
        first = regexprep(first, "[^ACGTU]", "");
        first = replace(first, "U", "T");

        if strlength(first) == 6
            tok(end+1,1) = first; %#ok<AGROW>
        end
    end

    hex = unique(tok, "stable");

    if isempty(hex)
        error("No valid hexamers parsed from: %s (check file format)", fp);
    end
end

function feat = computeMotifFeatures_ESE_ESS(seqs, motifs)
    n = numel(seqs);
    ESE_count = zeros(n,1);
    ESS_count = zeros(n,1);
    ESE_density = zeros(n,1);
    ESS_density = zeros(n,1);

    len = strlength(seqs);

    % Set true if you want "unique-position" style counting (optional)
    uniqueMode = false;

    for i = 1:n
        s = seqs(i);
        if len(i) == 0, continue; end

        if uniqueMode
            ESE_count(i) = countHitsUniquePositions(s, motifs.ESE);
            ESS_count(i) = countHitsUniquePositions(s, motifs.ESS);
        else
            ESE_count(i) = countHitsTotalOccurrences(s, motifs.ESE);
            ESS_count(i) = countHitsTotalOccurrences(s, motifs.ESS);
        end

        ESE_density(i) = 100 * ESE_count(i) / double(len(i));
        ESS_density(i) = 100 * ESS_count(i) / double(len(i));
    end

    ESE_minus_ESS    = ESE_density - ESS_density;
    ESE_to_ESS_ratio = (ESE_density + eps) ./ (ESS_density + eps);

    feat = table(len, ESE_count, ESS_count, ESE_density, ESS_density, ESE_minus_ESS, ESE_to_ESS_ratio);
end

function k = countHitsTotalOccurrences(seq, motifList)
    k = 0;
    for j = 1:numel(motifList)
        k = k + count(seq, motifList(j));
    end
end

function k = countHitsUniquePositions(seq, motifList)
    L = 6; % hexamer length
    n = strlength(seq);
    if n < L, k = 0; return; end

    hitPos = false(1, double(n));
    s = char(seq);

    for j = 1:numel(motifList)
        m = char(motifList(j));
        idx = strfind(s, m);
        for t = 1:numel(idx)
            hitPos(idx(t):idx(t)+L-1) = true;
        end
    end

    k = sum(hitPos);
end

function [drug,gene,cls] = parseContrastFilename(fp)
    [~,name,~] = fileparts(fp);
    parts = split(string(name), "_");
    drug = parts(1);
    gene = parts(2);
    cls  = parts(3);
end

function saveEnrichmentSummaries(outDir)
    combos = ["RG_SMN2","RG_ELP1","NVS_SMN2","NVS_ELP1"];

    for c = 1:numel(combos)
        cg = combos(c);
        drug = extractBefore(cg,"_");
        gene = extractAfter(cg,"_");

        respFile = fullfile(outDir, sprintf("%s_%s_responsive_sequence_functional.xlsx", drug, gene));
        nonrFile = fullfile(outDir, sprintf("%s_%s_nonresponsive_sequence_functional.xlsx", drug, gene));

        if ~exist(respFile,"file") || ~exist(nonrFile,"file")
            continue;
        end

        R = readtable(respFile);
        N = readtable(nonrFile);

        vars = ["ESE_count","ESS_count","ESE_density","ESS_density","ESE_minus_ESS","ESE_to_ESS_ratio", ...
                "sf2_hits","sf2_max","sc35_hits","sc35_max","srp40_hits","srp40_max","srp55_hits","srp55_max", ...
                "donor_mismatch_count","donor10_GC","exon7_GC","upIntron200_GC","downIntron200_GC"];

        vars = vars(ismember(vars, string(R.Properties.VariableNames)) & ismember(vars, string(N.Properties.VariableNames)));

        Summary = table(vars', 'VariableNames', {'Feature'});
        Summary.mean_responsive = nan(height(Summary),1);
        Summary.mean_nonresponsive = nan(height(Summary),1);
        Summary.log2FC_mean_Resp_vs_Nonresp = nan(height(Summary),1);
        Summary.p_value_ranksum = nan(height(Summary),1);

        for j=1:numel(vars)
            v = vars(j);
            r = R.(v); n = N.(v);
            r = r(~isnan(r)); n = n(~isnan(n));

            Summary.mean_responsive(j) = mean(r);
            Summary.mean_nonresponsive(j) = mean(n);
            Summary.log2FC_mean_Resp_vs_Nonresp(j) = log2(mean(r)+eps) - log2(mean(n)+eps);

            if ~isempty(r) && ~isempty(n)
                Summary.p_value_ranksum(j) = ranksum(r, n);
            end
        end

        outName = sprintf("%s_%s_ENRICHMENT_summary.xlsx", drug, gene);
        outPath = fullfile(outDir, outName);
        writetable(Summary, outPath);
        fprintf("Saved enrichment summary: %s\n", outPath);
    end
end

function [hitCounts, maxScores] = scanPWM(seqs, pwm, threshold)
    % Accept string OR char OR cellstr safely
    if isstring(seqs)
        seqs = cellstr(seqs);
    elseif ischar(seqs)
        seqs = cellstr(seqs); % single char row to cellstr
    end

    L = size(pwm, 2);
    hitCounts = zeros(numel(seqs),1);
    maxScores = nan(numel(seqs),1);

    for i = 1:numel(seqs)
        s = upper(seqs{i});
        s = strrep(s,'U','T');
        s = regexprep(s,'[^ACGT]','');
        n = length(s);

        if n < L
            hitCounts(i) = 0;
            maxScores(i) = NaN;
            continue;
        end

        scores = -inf(1, n-L+1);

        for k = 1:(n-L+1)
            w = s(k:k+L-1);
            sc = 0;
            ok = true;

            for p = 1:L
                switch w(p)
                    case 'A', r = 1;
                    case 'C', r = 2;
                    case 'G', r = 3;
                    case 'T', r = 4;
                    otherwise
                        ok = false;
                        break;
                end
                sc = sc + pwm(r,p);
            end

            if ok
                scores(k) = sc;
            end
        end

        maxScores(i) = max(scores);
        hitCounts(i) = sum(scores >= threshold);
    end
end