function [naiveoverallentropy,repoverallentropy]=calc_overallentropy(prelesionseq,postlesionseq,prelesionchunkseq,postlesionchunkseq)
% this function calculates overall transition entropy for sequences before
% and after the chi sq analysis. I have made it for the MMAN paper so that I
% can easily compare entropies before and after lesion
% inputs: prelesionseq, postlesionseq are naive sequences, usually with
% repeats removed (because for chi sq analysis you have to remove repeats
% and also that's what they do in the katahira 2013 paper. repeats = syl
% repeating > 2 times).
% prelesionchunkseq and postlesionchunkseq are outputs from the
% ak_chunkextraction function which gives you the seq with states replaced,
% but NOT chunks replaced. repeats are removed before giving the input to
% this function
% transition entropy is for all syllables, not just the branch points

% calculate overall entropy
% first use un-replaced, naive sequences
syls=unique(prelesionseq);
for ij=1:length(syls)
    labelspre{ij}=syls(ij);
end
[~,~,rawmtx]=calctransitionprob_fromsequence(labelspre,prelesionseq,1);
[~, preoverallte] = transentropy(rawmtx);
syls=unique(postlesionseq);
for ij=1:length(syls)
    labelspost{ij}=syls(ij);
end
[~,~,rawmtx]=calctransitionprob_fromsequence(labelspost,postlesionseq,1);
[~, postoverallte] = transentropy(rawmtx);
naiveoverallentropy=[preoverallte,postoverallte];

% now use the chi-sq analysis, replaced sequences
syls=unique(prelesionchunkseq);
for ij=1:length(syls)
    labelspre{ij}=syls(ij);
end
[~,~,rawmtx]=calctransitionprob_fromsequence(labelspre,prelesionchunkseq,1);
[~, r_preoverallte] = transentropy(rawmtx);

syls=unique(postlesionchunkseq);
for ij=1:length(syls)
    labelspost{ij}=syls(ij);
end
[~,~,rawmtx]=calctransitionprob_fromsequence(labelspost,postlesionchunkseq,1);
[~, r_postoverallte] = transentropy(rawmtx);
repoverallentropy=[r_preoverallte,r_postoverallte];
end