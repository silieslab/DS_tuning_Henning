function [mm, msglen] = mm_TryAddClass(mm, data)

if(mm.opts.display)
    fprintf('             Attempting to add and re-assign...\n');
end  

K = mm.nClasses + 1;

%% Add a new class
mm.class{K} = mm_CreateClass(mm.ModelTypes);
mm.nClasses = K;

% run kmeans 
warning('off', 'stats:kmeans:FailedToConverge');    % supress kmeans++ warnings
warning('off', 'stats:kmeans:MissingDataRemoved');
[~,~,~,Dist] = kmeans(data, K);
warning('on', 'stats:kmeans:FailedToConverge');   
warning('on', 'stats:kmeans:MissingDataRemoved');

R  = bsxfun(@rdivide, Dist, sum(Dist,2));
ix = sum(isnan(Dist),2) > 0;
if(any(ix))                
    t       = rand(sum(ix), K);
    R(ix,:) = bsxfun(@rdivide, t, sum(t,2));
end

mm.r  = R;
mm.Nk = sum(mm.r,1)';
mm.a  = mm.Nk / sum(mm.Nk);

% mm.a = zeros(K,1);
% mm.r = ones(size(data,1),K) * 0.05;
% for k = 1:K
%     ix = (IDX == k);
% 
%     % proportions
%     mm.a(k) = mean(ix);
%     mm.r(ix,k) = 0.95;
% end
% 
% 
% mm.r  = bsxfun(@rdivide, mm.r, sum(mm.r,2));
% mm.Nk = sum(mm.r,1)';     

%% Re-estimate the model parameters for the two classes using the new memberships
if(min(mm.Nk) >= mm.MinMembers)
    mm = mm_EstimateTheta(mm, data, 1:mm.nClasses);
end

%% Run the EM procedure
mm = mm_EM(mm, data);
msglen = mm.msglen;

if(mm.opts.display)
    fprintf('             Best      [ADD]: msglen = %10.2f nits\n', msglen);
end  

end