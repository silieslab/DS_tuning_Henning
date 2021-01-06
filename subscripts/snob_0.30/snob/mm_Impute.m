function data = mm_Impute(mm, data)

%% Determine the order of imputation - must do linear regression columns last
imputeOrder = [];
doLast = [];
for i = 1:mm.nModelTypes
    if( strcmp(mm.ModelTypes{1}.type,'linreg') )
        doLast = [doLast, i];
    end
    imputeOrder = [imputeOrder, i];
end
imputeOrder = [imputeOrder, doLast];

%% Do the imputation
for j = imputeOrder
    cols = mm.ModelTypes{j}.Ivar;    % cols to impute
    if(sum(sum(isnan(data(:,cols)))) > 0)
        data(:,cols) = mm_ImputeModel(mm, data(:,cols), j ,data);
    end
end

end

