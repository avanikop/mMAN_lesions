function [ent,histdepall]=transent_prevspost(prelesionseq,bpprelesion,postlesionseq,bppostlesion,lenv)
    % in prelesion seq find all of these branchpoints
    % calculate all possible transitions at these branchpoints and calculate
    % the raw numbers
    % using the raw numbers calculate transition entropy per branchpoint;
    % if you're giving regularexp anywhere, specify all lengths of
    % expressions in the lenv vector
    % old code: 
    % probmtx = rawmtx./repmat(sum(rawmtx,2),1,size(rawmtx,2));
    % rawmtx(probmtx<0.01) = 0;
    % sylcounts = sum(rawmtx,2);
    % 
    % probmtx = rawmtx./repmat(sylcounts,1,size(rawmtx,2));
    % 
    % 
    % te = nansum(-probmtx.*log2(probmtx),2)
    % I am also including history dependence here; should probably rename
    % the function later
    ent=zeros(length(bpprelesion),2);
    histdepall=zeros(length(bpprelesion),2);
    for ij=1:length(bpprelesion)
        if exist('lenv','var')
            ent(ij,1)=transitionentropy_bp(prelesionseq,bpprelesion{ij},lenv(ij));
            if ent(ij,1)==0
                histdepall(ij,1)=0;
            else
                histdepall(ij,1)=historydependence_bp(prelesionseq,bpprelesion{ij},lenv(ij));
            end
        else
            ent(ij,1)=transitionentropy_bp(prelesionseq,bpprelesion{ij});
            if ent(ij,1)==0
                histdepall(ij,1)=0;
            else
                histdepall(ij,1)=historydependence_bp(prelesionseq,bpprelesion{ij});
            end
        end
        if exist('postlesionseq','var') && exist('bppostlesion','var')
            if exist('lenv','var')
                ent(ij,2)=transitionentropy_bp(postlesionseq,bppostlesion{ij},lenv(ij));
                if ent(ij,2)==0
                    histdepall(ij,2)=0;
                else
                    histdepall(ij,2)=historydependence_bp(postlesionseq,bppostlesion{ij},lenv(ij));
                end
            else
                ent(ij,2)=transitionentropy_bp(postlesionseq,bppostlesion{ij});
                if ent(ij,2)==0
                    histdepall(ij,2)=0;
                else
                    histdepall(ij,2)=historydependence_bp(postlesionseq,bppostlesion{ij});
                end
            end
        end
    end

end 

