% Based on the function by Sveen Haase 
% https://au.mathworks.com/matlabcentral/fileexchange/23077-weighted-median
function wMed = medianw(x, w)

% normalize the weights
w = w / sum(w);

% sort the vectors
ASort = sortrows([x w],1);
dSort = ASort(:,1);
wSort = ASort(:,2);
sumVec = cumsum(wSort);    % vector for cumulative sums of the weights

% weighted median
wMed = dSort(find(sumVec>=0.5,1));

end