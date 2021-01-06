%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cn = removecells(c, index)

c(index) = [];
cn = cell(1, length(c)-sum(index));
j = 1;
for i = 1:length(c)
    if (~isempty(c{i}))
        cn{j} = c{i};
        j = j+1;
    end
end

end