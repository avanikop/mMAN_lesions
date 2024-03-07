function [con]=chunkconsistency(seq,expr,len)
%specify length of expr if you're using something like '[lB].'
% find consistency per chunks 14.07.22 Avani
if ~exist('len','var')
    len=length(expr);
end
point=(regexp(seq,expr));
%lena version
list2 = [];
for i = 1:length(point)
    list2 = [list2; seq(point(i):(point(i)+len-1))];
end
list=list2;
%what should Y char be? It should be simply LESS in double value than all
%unique syls in the sequence.
unq=unique(seq);
doubs=double(unq);
mindoubs=min(doubs);
chardoub=mindoubs-1;
Ychar=char(chardoub);

Yarray=[];
for ij=1:length(point')
    Yarray(ij)=chardoub;
end

newlist=[Yarray' list];
longstr='';
for ik=1:length(newlist)
    longstr=[longstr newlist(ik,:)];
end

%convert to cell array for labels
labels={};
unq=(unique(longstr));
for ix=1:length(unq)
    labels{ix}=unq(ix);
end
[a,divprob,TPnew]=calctransitionprob_fromsequence(labels,longstr,0);
% remove first row and first column
TPnew2= TPnew(2:end,2:end);
con = sum(max(TPnew2,[],2))./sum(sum(TPnew2));