%% Example - Mixture logistic regression - circular decision boundary
clear;

% Seed the random number generator
rng(7);

%% Generate data
n = 200;
x = rand(n, 2)*2 - 1;   % x \in [-1,1]
decr = .75;             % decision boundary radius
y = 1*(sqrt(x(:,1).^2 + x(:,2).^2) > decr);

%% Add noise
ix = randperm(n);
noise = round(n / 10);
y(ix(1:noise)) = abs(1 - y(ix(1:noise)));   

mm = snob([x, y], {'norm',1:2,'logreg',[3,1,2]}, 'k', 3);

% Print out the details of the model
mm_Summary(mm);

% Plot decision boundary
plot_db(mm,x,y);