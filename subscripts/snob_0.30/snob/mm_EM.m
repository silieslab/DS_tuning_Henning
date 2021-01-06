function mm = mm_EM(mm, data)

%% do search
done = false;
iter = 1;
L = zeros(mm.opts.emmaxiter, 1);

while not(done)
    
    %% Collapse mixtures if required ...
    if(min(mm.Nk) < mm.MinMembers)
        if(mm.opts.display)
            fprintf('             Removing %d class(es) due to insufficient membership...\n', sum(mm.Nk < mm.MinMembers));
        end                 
        mm = mm_Collapse(mm, mm.Nk < mm.MinMembers);
    end
    
    %% Estimate class proportions
    [mm.a, mm.r, mm.Nk, p] = mm_EstimateMixingWeights(mm, data);
        
    if(min(mm.Nk) >= mm.MinMembers)
        %% Estimate models parameters
        mm = mm_EstimateTheta(mm, data, 1:mm.nClasses);

        % Negative log-likelihood
        L(iter) = -sum(log(sum(p,2)));
        
        %% Check exit conditions
        cond = (iter > 5) && mean( abs(diff(L(iter-4:iter))) ) < 1e-2;
        if iter >= mm.opts.emmaxiter || cond
            done = true;
        end                

        % next iteration
        iter = iter + 1;
    end
    
end

%% Get the final message length of the model
mm = mm_MsgLen(mm, data);

end
