function [te]=transitionentropy_bp(seq,bp,len)
    % function for transition entropy for a given seq and branchpoint pair.
    % This function will then be used in ak_transent_postl.
    % can also choose to make this a separate function file
    if ~exist('len','var')
        len = length(bp);
    end
    try
        nextsylstr=seq(regexp(seq,bp)+len); %so that you can put a regexp
    catch
        try
            seq=seq(1:end-len);
            nextsylstr=seq(regexp(seq,bp)+len);
        catch
            nextsylstr=[];
        end
    end
    %so that you can add bp with
    % context of prev syl eg'ab' instead of 'b'
    if ~(length(nextsylstr)==0)
        alltrans=unique(nextsylstr);
        counts=zeros(1,length(alltrans));
        for i=1:length(alltrans)
            str=[bp,alltrans(i)];
            counts(i)=length(regexp(seq,str));
        end
        sylcounts=sum(counts);
        probmtx = counts./repmat(sylcounts,1,size(counts,2));
        te = sum(-probmtx.*log2(probmtx),2);
    else te=0;
    end
end