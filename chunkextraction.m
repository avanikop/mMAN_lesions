function[paths]=chunkextraction(g,thresh)
% find continuous path of 80% nodes
% seprate nodes as one in one out, multiple in multiple out,
% one in multiple out, multiple in one out
oimo=[];
mioo=[];
mimo=[];
oioo=[];
for i=1:g.numnodes
    [~,nIN]=inedges(g,i);
    [~,nOUT]=outedges(g,i);
    if length(nIN)==1 & length(nOUT)>1
        oimo=[oimo,i];
    elseif length(nIN)==1 & length(nOUT)==1
        oioo=[oioo,i];
    elseif length(nIN)>1 & length(nOUT)==1
        mioo=[mioo,i];
    else mimo=[mimo,i];
    end
end
nexplored=[];
paths={};
% use nodes in mioo (mult in one out) as starter nodes
for ij=1:length(mioo)
    v=mioo(ij); %also consider starting at 1
    [eid,w]=outedges(g,v);
    chkpath=[v];
    nexplored=[nexplored,v]; 
    u=w;
    if ismember(w,oioo)
        while g.Edges(eid,:).Weight>=thresh && ismember(w,oioo) %only oioo paths %TODO: oimo paths
            chkpath=[chkpath,w];
            [eidu,u]=outedges(g,w);
            nexplored=[nexplored,w];
            w=u;
            eid=eidu;     
        end
        nexplored=[nexplored,u]; %its okay if some us are double
        if g.Edges(eid,:).Weight>=thresh && ismember(u,oimo)
            chkpath=[chkpath,u];
        end
        paths=[paths,chkpath];
    elseif g.Edges(eid,:).Weight>=thresh && ismember(w,oimo)
        nexplored=[nexplored,w];
        chkpath=[chkpath,u];
        paths=[paths,chkpath];
    end
end
%% for chunks starting from v=Y
newset=setdiff(oioo,nexplored); %unexplored oioo nodes
name_start=findnode(g,'Y'); %%here you find the num of node Y
% I think I should force name_start to be considered here
newset=[newset,name_start];
newset=unique(newset); %in case name_start was already in the set

for ik=1:length(newset) %are sorted by seq anyway
    v=newset(ik);
    [~,k]=inedges(g,v);
    if v==name_start %%here
       [eid,w]=outedges(g,v);
       chkpath=[v];
       nexplored=[nexplored,v];
       u=w;
        if ismember(w,oioo)
            while g.Edges(eid,:).Weight>=thresh && ismember(w,oioo) %only oioo paths %TODO: oimo paths
                chkpath=[chkpath,w];
                [eidu,u]=outedges(g,w);
                nexplored=[nexplored,w];
                w=u;
                eid=eidu;     
            end
            nexplored=[nexplored,u]; %its okay if some us are double
            if g.Edges(eid,:).Weight>=thresh & ismember(u,oimo)
            chkpath=[chkpath,u];
            end
            paths=[paths,chkpath];
        elseif g.Edges(eid,:).Weight>=thresh & ismember(w,oimo)
            nexplored=[nexplored,w];
            chkpath=[chkpath,u];
            paths=[paths,chkpath];
        end
    end
    if v~=name_start && ~ismember(k,oioo) & ~ismember(v,nexplored) %%here
        [eid,w]=outedges(g,v);
        chkpath=[v];
        nexplored=[nexplored,v];
        u=w;
        if ismember(w,oioo)
            while g.Edges(eid,:).Weight>=thresh & ismember(w,oioo) %only oioo paths %TODO: oimo paths
                chkpath=[chkpath,w];
                [eidu,u]=outedges(g,w);
                nexplored=[nexplored,w];
                w=u;
                eid=eidu;     
            end
            nexplored=[nexplored,u]; %its okay if some us are double
            if g.Edges(eid,:).Weight>=thresh & ismember(u,oimo)
            chkpath=[chkpath,u];
            end
            paths=[paths,chkpath];
        elseif g.Edges(eid,:).Weight>=thresh & ismember(w,oimo)
            nexplored=[nexplored,w];
            chkpath=[chkpath,u];
            paths=[paths,chkpath];
        end
    end
end