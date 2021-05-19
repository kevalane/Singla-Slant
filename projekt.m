%% Simuleringsbaserat test
% Number of simulations
numberOfSims = 1e4;

% Fill a vector y with zeros of size numberOfSims
y = zeros(numberOfSims, 1);

% For loop to fill vector y
for i=1:numberOfSims
    % Simulate 200 independent tries
    z = binornd(1, 1/2, 2000, 1);
    
    % Append result to y-vector
    y(i) = find_sequence(z);
end

% Gives the mean of the distribution
mean(y, 'all')

% Plot result 
figure(1)
histogram(y, 'norm', 'prob', 'binmethod', 'int')

% Settings for the histogram
grid on;
xlabel('Svitlängd (st)');
ylabel('Sannolikhet (%)');
title('Histogram med de längsta sviterna per försök');

% And the cumulative distr
figure(2)
histogram(y, 'norm', 'cdf', 'binmethod', 'int', 'displaystyle', 'stair')

% Settings for the histogram
grid on;
xlabel('Svitlängd (st)');
ylabel('Sannolikhet (%)');
title('Kumulativ distribution av svitlängdernas utfall');


%% Regression
% Number of simulations
numberOfSims = 1e2;

% X vector of linspace over interesting x-values
% x = number of coin tosses in each simulation
%x = round( linspace(1e1, 1e3, 25) )';

% Using log space instead
x = round(logspace(1,4,25))';

% Fill y-vector with zeros of length x
y = zeros(length(x), 1);

% Loop over every single x-value top obtain y
for i=1:length(x)
    
    % Fill n vector with size number of sims
    n = zeros(numberOfSims,1);
    
    % Nested for loop to try each x-value numberOfSims times
    for j=1:numberOfSims
        % Get binomial where x(i) is from linspace
        z = binornd(1, 1/2, x(i), 1);
        % Fill n vector with longest suite
        n(j) = find_sequence(z);
    end
    % Fill corresponding y value (y(x)) with 
    % mean length of suite for each sim
    y(i) = mean(n);
end
% New figure
figure(3)

% Transform x-axis
x = log(x);

% Plot x to y
plot(x,y, 'o')
grid on;
xlabel('log(Antal slantsinglingar) (log(st))');
ylabel('Medelvärde svitlängd (st)');
title('Regressionslinje transformerad');

%% Parameter calculations
% Calculate Sxx 
Sxx = sum((x-mean(x)).^2);
% Calculate Sxy 
Sxy = sum((x-mean(x)).*(y-mean(y)));
% Calculate Syy 
Syy = sum((y-mean(y)).^2);

% Calculate Beta*
Beta = Sxy/Sxx;
% Calculate Alpha
ybar = mean(y);
xbar = mean(x);
Alfa = ybar - Beta*xbar;

% Calculate Q0
Q0 = Syy-(((Sxy)^2)/Sxx);
% Calculate deviation
s = sqrt(Q0/(length(n)-2));

% Calc deviation on B*
deviationBeta = s/sqrt(Sxx);
% Calc deviation on a*
deviationAlfa = s*sqrt(1/length(n) + (mean(x)^2)/Sxx);

% Deviation for prediction interval
deviationYofX = s*sqrt(1+(1/length(n)) + (((x - mean(x)).^2)./Sxx));

% Deviation for confidence interval of regression
deviationMuOfX = s*sqrt((1/length(n)) + (((x - mean(x)).^2)./Sxx));

% Deviation of calibration interval
deviationX0 = ((s/Beta)*sqrt(1 + 1/length(n) + (((y - mean(y)).^2)/((Beta^2)*Sxx))));

%% Plotting the regression line
% Getting the regression line
Y = Alfa + Beta*x;

% Plotting the regression line
hold on
plot(x, Y)

% Plot the mu lines 
hold on
muPlus = Y + 1.98.*deviationMuOfX;
muNeg = Y - 1.98.*deviationMuOfX;
plot(x, muPlus, x, muNeg)
title('Regressionslinje med konfidensintervall')

%% Calibration interval for 26 suite
% We want length 26 of suite
y0 = 26;
% Determine x0
x0 = (y0 - Alfa)/Beta;
% Deviation calculation
dX0 = ((s/Beta)*sqrt(1 + 1/length(n) + ((y0 - mean(y))^2)/((Beta^2)*Sxx)));
% Using t0.025(100)
err = 1.98*dX0;
% Calibration interval
Int = [x0 - err, x0 + err];
% Transform interval
Int = exp(Int)

