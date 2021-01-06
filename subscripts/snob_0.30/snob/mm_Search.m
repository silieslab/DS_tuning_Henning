function mm_best = mm_Search(mm, data)

%% Keep track of best model
mm_best = mm;
mm_best.msglen = inf;

%% possible outcomes of a split/merge attempt
outcomes = {'em','re-assign','split','merge','add'};

%% do search
done = false;
iter = 1;
msglen = zeros(mm.opts.maxiter,1);
prev_msglen = 0;
while not(done)

        %% Estimate model parameters using EM, given the model structure
        mm = mm_EM(mm, data);
        
        %% Try structural changes?
        if(~mm.opts.fixedstructure)
            [mm_s, msglen_s] = mm_TrySplits(mm, data);      % split           
            [mm_c, msglen_c] = mm_TryCombine(mm, data);     % combine
            [mm_a, msglen_a] = mm_TryAddClass(mm, data);    % add new class and re-assign
        else
            mm_c={}; mm_s={}; mm_a={}; msglen_s=inf; msglen_c=inf; msglen_a=inf;
        end        
        
        %% Attempt to re-assign items to classes
        if(mm.nClasses > 1)
            [mm_r, msglen_r] = mm_Reassign(mm, data);        
        else
            mm_r = {}; msglen_r = inf;
        end

        %% Pick model stochastically based on posterior probability
        mm_cand = {mm; mm_r; mm_s; mm_c; mm_a};
        msglen_cand = [mm.msglen; msglen_r; msglen_s; msglen_c; msglen_a];
        
        if (~mm.opts.greedy)
            % Get induced probability distribution over all candidates
            mm_prob = msglen_cand - min(msglen_cand);
            mm_prob = exp(-mm_prob)/sum(exp(-mm_prob));

            % Sample a model
            I = mnrnd(1, mm_prob) == 1;
            mm = mm_cand{I};

        %% Otherwise greedy search (pick model with smallest message length)
        else
            [~,I] = min(msglen_cand);
            mm = mm_cand{I};
        end
                
        %% Is it the best model so far?
        if (mm.msglen < mm_best.msglen)
            mm_best = mm;   % keep track of the best model
        end
        
        %% Record current best msglen
        msglen(iter) = mm_best.msglen;
        
        if(mm.opts.display)
            fprintf('             OUTCOME: %s [%+10.2f]\n', outcomes{I}, msglen(iter) - prev_msglen);
        end             
        prev_msglen = msglen(iter);
        
        %% Check exit conditions
        cond = (iter > 5) && mean( abs(diff(msglen(iter-4:iter))) ) < 1e-2;
        if iter >= mm.opts.maxiter || cond
            done = true;
        end      
        
        %% Print current state
        if(mm.opts.display)
            fprintf('[Iter %6d: #classes = %2d, msglen = %10.2f nits, L = %10.2f, cost = %5.2f]\n', iter, mm_best.nClasses, mm_best.msglen, mm_best.L, mm_best.Atheta + mm_best.constant + mm_best.Ak + mm_best.Aa);
        end        
        
        % next iteration
        iter = iter + 1;
    
end

end