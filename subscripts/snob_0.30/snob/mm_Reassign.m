function [mm_r, msglen] = mm_Reassign(mm, data)

mm_r = mm;

%% randomly re-assign all items
r = mm.r;
n = size(data,1);
mm_r.r = r(randperm(n),:);
mm_r.Nk = sum(r, 1)';         

% initial estimates of theta based on the re-assignment
mm_r = mm_EstimateTheta(mm_r, data, 1:mm_r.nClasses);   

% optimise the current model
mm_r = mm_EM(mm_r, data);
msglen = mm_r.msglen;

if(mm.opts.display)
    fprintf('               Re-assignment: msglen = %10.2f nits\n', mm_r.msglen);
end   

end