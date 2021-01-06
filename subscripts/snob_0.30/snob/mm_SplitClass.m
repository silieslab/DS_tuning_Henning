function mm = mm_SplitClass(mm, data, i)

n = size(data,1);
K = mm.nClasses;

%% Add a new class
mm.class{K+1} = mm_CreateClass(mm.ModelTypes);
mm.nClasses = K+1;

%% Randomly divide the memberships of the datapoints between the classes
mix = rand(n,1);
r = mm.r(:,i);
mm.r(:,i) = r.*mix;
mm.r(:,mm.nClasses) = r.*(1-mix);

mm.Nk = sum(mm.r,1)';
mm.a = (mm.Nk+1/2)./(n+mm.nClasses/2);

%% Re-estimate the model parameters for the two classes using the new memberships
if(min(mm.Nk) >= mm.MinMembers)
    mm = mm_EstimateTheta(mm, data, [i, mm.nClasses]);
end

%% Run the EM procedure
mm = mm_EM(mm, data);

end