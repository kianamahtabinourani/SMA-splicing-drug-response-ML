%% 00 Extract Smn2 Exon7
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


inputFasta  = fullfile(projectRoot, "data", "raw", "SMN2_ensembl.fasta");
outputFasta = fullfile(projectRoot, "data", "intermediate", "SMN2_exon7_200bp.fasta");

% ✅ Prevent appending: delete old output if it exists
if exist(outputFasta, 'file')
    delete(outputFasta);
end

seqs = fastaread(inputFasta);

% ---- Find exon 7 record by Ensembl exon ID ----
targetExonID = 'ENSE00003698693';
exon7Seq = '';

for i = 1:length(seqs)
    if contains(seqs(i).Header, targetExonID)
        exon7Seq = upper(seqs(i).Sequence);
        break;
    end
end

if isempty(exon7Seq)
    error('Exon 7 sequence not found by exon ID.');
end

fprintf('Exon 7 length: %d bp\n', length(exon7Seq));

% ---- Find chromosome sequence record ----
chrSeq = '';
for i = 1:length(seqs)
    if startsWith(seqs(i).Header, '5 dna:chromosome')
        chrSeq = upper(seqs(i).Sequence);
        break;
    end
end

if isempty(chrSeq)
    error('Chromosome sequence record not found (">5 dna:chromosome ...").');
end

% ---- Locate exon 7 inside chromosome sequence ----
pos = strfind(chrSeq, exon7Seq);

if isempty(pos)
    error('Could not locate exon 7 sequence within chromosome sequence.');
end

startPos = pos(1);
endPos   = startPos + length(exon7Seq) - 1;

% ---- Extract ±200 bp ----
flank = 200;
leftStart = max(1, startPos - flank);
rightEnd  = min(length(chrSeq), endPos + flank);

finalSeq = chrSeq(leftStart:rightEnd);

fprintf('Final extracted length: %d bp\n', length(finalSeq));

% ---- Save ----
fastawrite(outputFasta, 'SMN2_exon7_plus_minus_200bp', finalSeq);

disp('✅ Done. Saved file:');
disp(outputFasta);