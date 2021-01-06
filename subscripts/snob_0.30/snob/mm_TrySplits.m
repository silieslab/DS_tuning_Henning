function [best_model, best_msglen] = mm_TrySplits(mm, data)

if(mm.opts.display)
    fprintf('             Attempting to split...\n');
end  

% number of mixtures; number of regressors
K = mm.nClasses;  

% keep track of best model
[~, II] = max(mm.r, [], 2);        

%% attempt to split each existing class 
mm_s = cell(K,1);
msglen_s = ones(K,1)*inf;
for i = 1:K
    
    ix = (II == i);
    
    % Do we have enough data in this class to split it?
    if( sum(ix) > 2 * mm.MinMembers )    
        mm_s{i} = mm_SplitClass(mm, data, i);
        msglen_s(i) = mm_s{i}.msglen;

        if(mm.opts.display)
            fprintf('             class %3d: msglen = %10.2f nits\n', i, msglen_s(i));        
        end
    end
end

%% Return the best split
if(all(isinf(msglen_s)))
    best_model = mm;
    I = 0;
else    
    [~,I] = min(msglen_s);
    best_model = mm_s{I};
end
best_msglen = best_model.msglen;

if(mm.opts.display)
    fprintf('             Best      [%3d]: msglen = %10.2f nits\n', I, best_msglen);
end     

end