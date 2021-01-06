function plot_db(mm, X, y)

%% Plot the points

%% Plot decision plane
x1range = linspace(min(X(:,1)),max(X(:,1)),1e3);
x2range = linspace(min(X(:,2)),max(X(:,2)),1e3);
[xx1, xx2] = meshgrid(x1range,x2range);
XGrid = [xx1(:) xx2(:)];

%% Compute predictions [assume logreg model is last]
n = size(XGrid,1);
ypred=zeros(n,mm.nClasses); 
for i = 1:mm.nClasses
    b0 = mm.class{i}.model{end}.theta(1);
    b  = mm.class{i}.model{end}.theta(2:end);
    ypred(:,i) = logsig(b0 + XGrid*b); 
end
r = mm_EstimateR(mm, [XGrid,nan(n,1)]);
ypred = sum( ypred .* r, 2) > 0.5;

%% Plot decision boundary and actual data 
gscatter(xx1(:), xx2(:), ypred, [0.4 0.4 0.4; 0.8 0.8 0.8],'.',10,'off');
hold;
gscatter(X(:,1),X(:,2),y,'br','ox',10,'off');
grid;
xlim([x1range(1),x1range(end)]);
ylim([x2range(1),x2range(end)]);
set(gcf,'color','white');
xlabel('x1','Fontsize',18);
ylabel('x2','Fontsize',18);
title('Decision boundary for mixture of logistic regression models', 'Fontsize', 16);

end
