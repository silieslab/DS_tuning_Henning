function [mm_combine, best_msglen] = mm_TryCombine(mm, data)

%% Only possible to combine if K > 1
K = mm.nClasses;
if (K > 1)
    if(mm.opts.display)
        fprintf('             Attempting to combine...\n');
    end            
    
    CombList = combnk(1:K, 2);
    nc = size(CombList ,1);
    
    %% If more combinations than allowed, randomly select a subset of them
    if (nc > mm.opts.MaxTryCombines)
        I = randperm(nc);
        CombList = CombList(I(1:mm.opts.MaxTryCombines),:);
        nc = mm.opts.MaxTryCombines;
    end
    
    %% Try all combinations
    mm_c = cell(nc,1);
    msglen_c = zeros(nc,1);
    for i = 1:nc       
        % Try combining
        mm_c{i} = mm_CombineClasses(mm, data, CombList(i,:));
        msglen_c(i) = mm_c{i}.msglen;
        
        if(mm.opts.display)
            fprintf('             class [%3d,%3d]: msglen = %10.2f nits\n', CombList(i,1), CombList(i,2), msglen_c(i));        
        end            
    end

    %% Return the best combination
    [~,I] = min(msglen_c);
    mm_combine = mm_c{I};
    best_msglen = mm_combine.msglen;
    
    if(mm.opts.display)
        fprintf('             Best  [%3d,%3d]: msglen = %10.2f nits\n', CombList(I,1), CombList(I,2), best_msglen);
    end         
    
%% Otherwise return empty lists    
else
    mm_combine = {};
    best_msglen = inf;
end



end