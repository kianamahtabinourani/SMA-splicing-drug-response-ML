%% 01 Extract Splice Regulatory Features
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% SMN2: Convert 5'SS variants into full splice-regulatory (sequence) features ONLY
% No PSI features are computed here.

clc; clear all; close all; 
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% 1) Paths (EDIT THESE)
smn2Csv    = fullfile(projectRoot, "data", "raw", "GSE221868_smn2_psi.csv");

% Your Ensembl FASTA: 200 nt upstream + exon7 + 200 nt downstream
exon7Fasta = fullfile(projectRoot, "SMN2_exon7_200bp.fasta"); 

outXlsx = fullfile(projectRoot, "Full splice regulatory dataset", "SMN2_splice_regulatory_features_ONLY.xlsx");

%% 2) Read SMN2 5'SS table
T = readtable(smn2Csv);

% ---- IMPORTANT ----
% If your first column name is not "ss", uncomment the next 2 lines:
% firstCol = T.Properties.VariableNames{1};
% T.Properties.VariableNames{1} = "ss";

ss = upper(string(T.ss));   % 10-nt 5' splice-site sequences

%% 3) Read exon7 ±200 FASTA (template context)
F = fastaread(exon7Fasta);
template = upper(string(F.Sequence));

% Expect 454 nt = 200 + 54 + 200 (for your file)
if strlength(template) ~= 454
    warning("Template length is %d (expected 454). Check exon boundaries below.", strlength(template));
end

UP   = extractBetween(template, 1,   200);
EXON = extractBetween(template, 201, 254); % SMN2 exon 7 (expected 54 nt)
DOWN = extractBetween(template, 255, strlength(template));

%% 4) Reference donor window from template: last 4 exon + first 6 intron = 10-mer
exonTail4_ref = extractBetween(EXON, strlength(EXON)-3, strlength(EXON));
downHead6_ref = extractBetween(DOWN, 1, 6);
donor10_ref   = exonTail4_ref + downHead6_ref;

%% 5) Build full sequence per 5'SS variant (substitute the 10-mer into the template)
n = height(T);

fullSeq        = strings(n,1);
exonSeq        = strings(n,1);
upIntron200    = strings(n,1);
downIntron200  = strings(n,1);
donor10        = strings(n,1);

for i = 1:n
    s = ss(i);
    if strlength(s) ~= 10
        error("Row %d: ss length is %d, expected 10.", i, strlength(s));
    end

    % Map 10-mer as: last 4 nt of exon + first 6 nt of downstream intron
    exonTail4 = extractBetween(s, 1, 4);
    downHead6 = extractBetween(s, 5, 10);

    % Replace last 4 nt of exon
    EXON_i = replaceBetween(EXON, strlength(EXON)-3, strlength(EXON), exonTail4);

    % Replace first 6 nt of downstream intron
    DOWN_i = replaceBetween(DOWN, 1, 6, downHead6);

    % Construct full context sequence
    fullSeq(i) = UP + EXON_i + DOWN_i;

    % Store parts
    donor10(i)       = exonTail4 + downHead6;
    exonSeq(i)       = EXON_i;
    upIntron200(i)   = UP;
    downIntron200(i) = DOWN_i;
end

%% 6) Compute splice-regulatory SEQUENCE features (simple + useful)
Out = table();
Out.ss = ss;

% donor sequence windows
Out.donor10_ref = repmat(donor10_ref, n, 1);
Out.donor10     = donor10;

% canonical donor check: intron starts with GT (RNA GU)
Out.has_GT = (extractBetween(donor10, 5, 6) == "GT");

% mismatches vs reference donor10
Out.donor_mismatch_count = zeros(n,1);
for i = 1:n
    Out.donor_mismatch_count(i) = sum(char(donor10(i)) ~= char(donor10_ref));
end

% GC% features
Out.donor10_GC     = gc_percent(donor10);
Out.exon7_GC       = gc_percent(exonSeq);
Out.upIntron200_GC = gc_percent(upIntron200);
Out.downIntron200_GC = gc_percent(downIntron200);

% simple pyrimidine fraction near acceptor region (last 30 nt of upstream intron)
% (proxy for polypyrimidine tract strength; 0..1)
upLast30 = extractBetween(upIntron200, 171, 200);
Out.upIntron_last30_pyFrac = pyrimidine_fraction(upLast30);

% Optional: save sequences too (comment these out if you want smaller output)
Out.fullSeq = fullSeq;
Out.exonSeq = exonSeq;
Out.upIntron200 = upIntron200;
Out.downIntron200 = downIntron200;

%% 7) Save
writetable(Out, outXlsx);
disp("Saved splice-regulatory features ONLY to: " + outXlsx);

%% ---- helper functions ----
function gc = gc_percent(seqStr)
    seqStr = upper(string(seqStr));
    gc = zeros(numel(seqStr),1);
    for k = 1:numel(seqStr)
        s = char(seqStr(k));
        if isempty(s), gc(k)=NaN; continue; end
        gc(k) = 100 * (sum(s=='G' | s=='C') / numel(s));
    end
end

function py = pyrimidine_fraction(seqStr)
    seqStr = upper(string(seqStr));
    py = zeros(numel(seqStr),1);
    for k = 1:numel(seqStr)
        s = char(seqStr(k));
        if isempty(s), py(k)=NaN; continue; end
        py(k) = (sum(s=='C' | s=='T') / numel(s));
    end
end