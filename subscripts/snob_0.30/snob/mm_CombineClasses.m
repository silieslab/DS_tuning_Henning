function mm = mm_CombineClasses(mm, data, pair)

%% Remove the class with the larger number
mm.class = removecells(mm.class, pair(2));

%% Remaining class gets all membership from the removed class
mm.r(:,pair(1)) = mm.r(:,pair(1)) + mm.r(:,pair(2));
mm.r(:,pair(2)) = [];
mm.nClasses = mm.nClasses - 1;

%% Update mixing proportions
mm.Nk = sum(mm.r,1)';
mm.a = (mm.Nk+1/2)./(size(data,1)+mm.nClasses/2);

%% Re-estimate the model parameters for the two classes using the new memberships
mm = mm_EstimateTheta(mm, data, pair(1));

%% Finally, run the EM procedure
mm = mm_EM(mm, data);

end