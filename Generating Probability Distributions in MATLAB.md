# Binomial Distribution
 The pmf of the binomial distribution is given by $$Bin(k|n, p) = \binom{n}{k}p^k(1 - p)^{n - k}$$
 To display it in MATLAB, we can use two approaches:
 - The first is to consider it as sum of $N$ iid bernoulli(p) random variables. Here we are using the property that if $(X_i)_{i = 1}^{N} \sim Ber(p)$ are iids then $\sum_{i = 1}^{N} X_i \sim Bin(n, p)$.
 - The second and more convenient approach is to use $binord$ in MATLAB.
We display the use of both approaches.

## The first approach

~~~matlab
% Set the parameter for the Bernoulli distribution
p = 0.1; % You can change this value as needed

% Number of repetitions
N = 1000;

% Number of Bernoulli random variables to generate each time
n = 20;

% Initialize an array to store the sums
sumsArray = zeros(1, N);

% Repeat the process N times
for i = 1:N
    % Generate n Bernoulli random variables
    berRandomVariables = rand(1, n) < p;
    
    % Add the generated random variables and store in the array
    sumsArray(i) = sum(berRandomVariables);
end

% Create a table to observe the frequency of elements from 1 to n
freqTable = tabulate(sumsArray);

% Display the table
disp('Frequency Table:');
disp('----------------');
disp('Value   Frequency');
disp('-----   ---------');
for i = 1:size(freqTable, 1)
    fprintf('%4d      %6.2f%%\n', freqTable(i, 1), freqTable(i, 3));
end

% Create a histogram
figure;
bar(freqTable(:, 1), freqTable(:, 3) / 100, 'BarWidth', 0.8);
xlabel('Observed Values');
ylabel('Probability');
title(['Bin(', num2str(n), ', ', num2str(p), ') Distribution']);
grid on;
~~~

## The Second Approach










# Poisson Distribution

$$Pois(k|\lambda) = e^{-\lambda} \frac{\lambda^k}{k!}$$
~~~matlab
% Parameters for the Poisson distribution
lambda = 3; % Average rate (parameter of the Poisson distribution)

% Generate a range of values for the random variable (e.g., from 0 to 10)
x_values = 0:100;

% Calculate the probabilities using the poisspdf function
probabilities = poisspdf(x_values, lambda);

% Plot the probabilities using a bar plot
bar(x_values, probabilities, 'b');
xlabel('X (Poisson Random Variable)');
ylabel('Probability');
title(['Pois(\lambda = ' num2str(lambda) ') Distribution']);

% Display the plot
grid on;
~~~



# Gaussian Distribution

~~~matlab
% Take user input for mean (u) and standard deviation (sigma)
u = input('Enter the mean (u): ');
sigma = input('Enter the standard deviation (sigma): ');

% Define the range of values
x_values = linspace(u - 4 * sigma, u + 4 * sigma, 1000);

% Calculate the PDF of the normal distribution
pdf_values = normpdf(x_values, u, sigma);

% Plot the PDF
plot(x_values, pdf_values, 'b-', 'LineWidth', 2);

% Label the axes and add a title
xlabel('X');
ylabel('Probability Density');
title(['PDF of N(' num2str(u) ', ' num2str(sigma) ')']);

% Display the plot
grid on;
~~~


# Student t Distribution

~~~matlab
% Values for degrees of freedom (nu)
degrees_of_freedom = [0.5, 1, 10, 100];

% Generate a range of values for the random variable (e.g., from -5 to 5)
t_values = linspace(-5, 5, 1000);

% Create a figure for the plot
figure;

% Plot the PDFs for different degrees of freedom on the same graph
hold on;
for i = 1:length(degrees_of_freedom)
    nu = degrees_of_freedom(i);
    pdf_values = tpdf(t_values, nu);
    
    plot(t_values, pdf_values, 'LineWidth', 2, 'DisplayName', ['\nu = ' num2str(nu)]);
end

% Label the axes and add a title
xlabel('t');
ylabel('Probability Density');
title('PDFs of Student''s t-distribution for Different Degrees of Freedom');
legend('show');
grid on;
hold off;
~~~

# The Cauchy Distribution

~~~matlab
% Parameters for the Cauchy distribution
x0 = 0; % Location parameter (median)
gamma = 1; % Scale parameter

% Generate a range of values for the random variable
x_values = linspace(-10, 10, 1000);

% Calculate the PDF of the Cauchy distribution manually
pdf_values = 1./(pi * gamma * (1 + ((x_values - x0)/gamma).^2));

% Plot the PDF
plot(x_values, pdf_values, 'b-', 'LineWidth', 2);

% Label the axes and add a title
xlabel('x');
ylabel('Probability Density');
title('PDF of Cauchy Distribution');

% Display the plot
grid on;
~~~



# The Gamma Distribution

~~~matlab
% Define parameters
a_values = [0.5, 1, 2];
b_values = [0.5, 1, 2];  % Different values of b

% Generate values for x
x = linspace(0, 10, 100);

% Plot the gamma PDF for different values of a and b
figure;
hold on;

for i = 1:length(a_values)
    for j = 1:length(b_values)
        a = a_values(i);
        b = b_values(j);
        pdf_gamma = gampdf(x, a, b);
        plot(x, pdf_gamma, 'LineWidth', 2, 'DisplayName', ['a = ' num2str(a) ', b = ' num2str(b)]);
    end
end

hold off;

% Add labels and legend
title('Gamma Distribution PDF for Different Values of a and b');
xlabel('x');
ylabel('Probability Density');
legend('Location', 'best');
grid on;
~~~


# The Beta Distribution

~~~matlab
% Parameters
a_values = [1, 0.1, 2, 2, 8];
b_values = [1, 0.1, 2, 3, 4];

% Plotting
figure;

hold on;

for i = 1:length(a_values)
    a = a_values(i);
    b = b_values(i);
    
    x = linspace(0, 1, 1000);
    y = betapdf(x, a, b);
    
    plot(x, y, 'LineWidth', 2, 'DisplayName', ['a = ' num2str(a) ', b = ' num2str(b)]);
end

hold off;

title('Beta PDFs for Different Parameter Sets');
xlabel('x');
ylabel('Probability Density');
legend('show');
grid on;

% Adjust y-axis limits
ylim([0, 5]); % Adjust the limits as needed
~~~


# The Pareto Distribution

~~~matlab
% Parameters
m_values = [0.01, 0.00, 1.00];
k_values = [0.10, 0.50, 1.00];

% Plotting
figure;

hold on;

for i = 1:length(m_values)
    m = m_values(i);
    k = k_values(i);
    
    x = linspace(m, 5, 1000);
    y = k * m.^k .* x.^(-k - 1) .* (x >= m);
    
    plot(x, y, 'LineWidth', 2, 'DisplayName', ['m = ' num2str(m) ', k = ' num2str(k)]);
end

hold off;

title('Pareto PDFs for Different Parameter Sets');
xlabel('x');
ylabel('Probability Density');
legend('show');
grid on;
~~~


# Plotting Contours for 2D Gaussian

~~~matlab
% Define the mean vector mu
mu = [0, 0];

% Define the covariance matrix Sigma
Sigma = [1, 0.5; 0.5, 2];

% Generate a grid of points
[x, y] = meshgrid(linspace(-5, 5, 100), linspace(-5, 5, 100));

% Reshape the grid into a column vector for efficient matrix operations
grid = [x(:), y(:)];

% Compute the Gaussian distribution values for each point in the grid
z = mvnpdf(grid, mu, Sigma);

% Reshape the output to match the shape of the input grid
z = reshape(z, size(x));

% Plot the 2D Gaussian distribution
figure;
contour(x, y, z, 'linewidth', 2);
title('Full Covariance matrix');
xlabel('X-axis');
ylabel('Y-axis');
~~~


~~~matlab
% Define the mean vector mu
mu = [0, 0];

% Define the covariance matrix Sigma
Sigma = [1, 0; 0, 2];

% Generate a grid of points
[x, y] = meshgrid(linspace(-5, 5, 100), linspace(-5, 5, 100));

% Reshape the grid into a column vector for efficient matrix operations
grid = [x(:), y(:)];

% Compute the Gaussian distribution values for each point in the grid
z = mvnpdf(grid, mu, Sigma);

% Reshape the output to match the shape of the input grid
z = reshape(z, size(x));

% Plot the 2D Gaussian distribution
figure;
contour(x, y, z, 'linewidth', 2);
title('Diagonal Covariance matrix');
xlabel('X-axis');
ylabel('Y-axis');
~~~


~~~matlab
% Define the mean vector mu
mu = [0, 0];

% Define the covariance matrix Sigma
Sigma = [1.5, 0; 0, 1.5];

% Generate a grid of points
[x, y] = meshgrid(linspace(-5, 5, 100), linspace(-5, 5, 100));

% Reshape the grid into a column vector for efficient matrix operations
grid = [x(:), y(:)];

% Compute the Gaussian distribution values for each point in the grid
z = mvnpdf(grid, mu, Sigma);

% Reshape the output to match the shape of the input grid
z = reshape(z, size(x));

% Plot the 2D Gaussian distribution
figure;
contour(x, y, z, 'linewidth', 2);
title('Spherical Covariance matrix');
xlabel('X-axis');
ylabel('Y-axis');
~~~

To plot the 2-D Gaussian

~~~matlab
% Define a grid of x and y values
[x, y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Compute the function values for each point in the grid
z = 1/sqrt(pi) * exp(-x.^2 - y.^2);

% Plot the surface
figure;
surf(x, y, z);

% Add labels and title
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('z = 2-D Gaussian');
~~~