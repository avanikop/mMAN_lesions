function [te, overallte] = transentropy(rawmtx)
% transition entropy per branch point
% and overall transition entropy, weighted by frequency of syllable
% as in Katahira 2013
% just te per bp should work with prob matrix as well

sylcounts = sum(rawmtx,2); % sum of all columns in a row

probmtx = rawmtx./repmat(sylcounts,1,size(rawmtx,2));


te = nansum(-probmtx.*log2(probmtx),2);

overallte = -sum(nansum(probmtx.*log2(probmtx),2).*(sylcounts./sum(sylcounts)));