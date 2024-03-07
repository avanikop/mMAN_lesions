function [outelts, outfreq] = uniquestring(seq,thresh)
% get unique elements in a string above a certain threshold (e.g. 2% of all
% elements)

[a,~,c]=unique(seq);
chist = histc(c,unique(c));
chist = chist./(sum(chist));
outelts = a(chist>thresh./100);
outfreq = chist(chist>thresh./100);