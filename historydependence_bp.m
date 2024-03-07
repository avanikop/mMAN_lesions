function [histdep]=historydependence_bp(seq,bp,len)
    % function for history dependence for a given seq and branchpoint pair
    % this function will be used in ak_transent_postl
    % could amke it into a separate function also
    % find the most probable transition 'ab'
    if ~exist('len','var')
        len = length(bp);
    end
    try
        nextsylstr=seq(regexp(seq,bp)+len); %so that you can put a regexp
    catch
        seq=seq(1:end-len);
        nextsylstr=seq(regexp(seq,bp)+len);
    end
    alltrans=unique(nextsylstr);
    counts=zeros(1,length(alltrans));
    for i=1:length(alltrans)
        str=[bp,alltrans(i)];
        counts(i)=length(regexp(seq,str));
    end
    sylcounts=sum(counts);
    probmtx = counts./repmat(sylcounts,1,size(counts,2));
    b=alltrans(probmtx==max(probmtx));
    findabadot=nextsylstr(strfind(nextsylstr(1:end-1),b)+1); %find all transitions after the dominant transition
    countabab=length(strfind(findabadot,b));%count how many times 'b' appears given that 'b' is the transitin at n-1
    probabab=countabab/length(findabadot); % prob of ab|ab = #ab|ab/#a-all|ab
    %for the second part, find all transitions that are NOT ab
    acpos=regexp(nextsylstr(1:end-1),['[^',b,']']); %all transitions that are NOT ab (here called ac)
    findacadot=nextsylstr(acpos+1); %find all transitions if ac is the transition previously
    countacab=length(strfind(findacadot,b)); % fing out how many of those are ab
    probacab=countacab/length(findacadot); % prob of ab|ac = #ab|ac/#a-all|ac

    % find prob of (ab at n|ab at n-1)
    % find prob of (ab at n|ac at n-1) where ac is any other transition
    % from a ie a-other
    % histdep=abs(p(ab at n|ac at n-1)-p(ab at n|ab at n-1))
    if isnan(probacab)
        probacab = 0;
    end
    if isnan(probabab)
        probabab = 0;
    end
    histdep=abs(probacab-probabab);

end