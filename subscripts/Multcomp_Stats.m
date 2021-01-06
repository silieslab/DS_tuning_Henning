function Multcomp_Stats(Data,Con)

% multiple comparison statistics

%Data and Con must have the same dimensions 

if size(Data,1)==length(Con)

    % first test data for normal distribution 
    Normal=nan(1,length(Con)); 
    for i=1:length(Con)
         Normal(i) = lillietest(Data(i,:));
    end 

    
    if sum(Normal)==0 % If all test groups are normally distributed 
        % Calc ANOVA
        [~,~,stats]=anova1(Data',[], 'off'); % It ignores NaNs!! (I tested with the group option, removing all NaNs and get the same statistics
        figure; [c,~,~,~] = multcompare(stats, 'CType', 'bonferroni','Display', 'on');     
        Test='ANOVA';
    else % If one test group is not normally distributed 
        % Calc Kruskalwallis for non normally distributed data
        [~,~,stats]=kruskalwallis(Data',[], 'off'); % It ignores NaNs!! (I tested with the group option, removing all NaNs and get the same statistics
        figure; [c,~,~,~] = multcompare(stats, 'CType', 'bonferroni','Display', 'on');
       Test='KKW';
    end 


     
          Pval=c(:,6);
        for iii=1:length(Con)-1
            N1=length(find(~isnan(Data(1,:))));
            N2=length(find(~isnan(Data(1+iii,:))));
            disp([Test, ': ', Con{1}, '(N=', num2str(N1),')- ' Con{1+iii}, '(N=', num2str(N2),') p-value: ', num2str(Pval(iii))])
        end 
   
     
        if size(Data,1)>3
            for iii=1:length(Con)-2
                N1=length(find(~isnan(Data(2,:))));
                N2=length(find(~isnan(Data(2+iii,:))));
                disp([Test, ': ', Con{2}, '(N=', num2str(N1),')- ' Con{2+iii}, '(N=', num2str(N2),') p-value: ', num2str(Pval(3+iii))])
            end 
        end 
        
        if size(Data,1)>4
            for iii=1:length(Con)-3
                N1=length(find(~isnan(Data(3,:))));
                N2=length(find(~isnan(Data(3+iii,:))));
                disp([Test, ': ', Con{3}, '(N=', num2str(N1),')- ' Con{3+iii}, '(N=', num2str(N2),') p-value: ', num2str(Pval(5+iii))])
            end 
        end 
else 
    disp('Dimensions of Condition and Data do not agree')
end 

end 