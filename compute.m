%% Labs in Numerical Modeling - Tutorial 1: Polynomial Approach, Group-5
close,clc,clear

%% 1. Interpolation case
% Given measurements for interpolation
xi_interp = [-5.0; 0.0; 10.0; 20.0; 25.0];
yi_interp = [4.8; 2.9; 3.5; 2.4; 3.6];

N = length(xi_interp) - 1; % the degree of the polynomial for interpolation case

%--- 1.a)
% Call the function coefficients.m to find the coefficients and the correlation matrix A
[c,A] = coefficients(xi_interp, yi_interp, N);

%--- 1.b)
% Evaluate the polynomial at x = 18
x_value = 18;
powers_x = x_value.^(0:N); % contains the powers of x value(18) from 0 to N
fx_18_value = dot(powers_x, c);

%% 2. Approximation case
% Add additional given measurements for approximation
additional_xi = [5.0; 15.0]; additional_yi = [3.1; 2.8];
xi_approx = [xi_interp; additional_xi];
yi_approx = [yi_interp; additional_yi];

%--- 2.a)
N = 4; % the degree of the polynomial for approximation case
% Calculate coefficients - LSE and A matrix
[c_approx, A_approx] = coefficients(xi_approx, yi_approx,N);

%--- 2.b)
% Evaluate the polynomial at the x value (18) using the dot product of powers_x_approx and c_approx
powers_x_approx = x_value.^(0:N); % contains the powers of x value(18) from 0 to N
fx_18_approx_value = dot(powers_x_approx, c_approx);

% Compare the result with the result from task 1.b)
disp('The value f(x) for x = 18 by Interpolation:'); disp(fx_18_value);
disp('The value f(x) for x = 18 by Approximation:'); disp(fx_18_approx_value);

%--- 2.c)
% Compare the calculated coefficients from the interpolation and the approximation case
disp('----For Comparison of Coefficients----'); disp('Interpolation Coefficients:'); disp(c);
disp('Approximation Coefficients:'); disp(c_approx);
% The difference between the two vectors of coefficients
coeff_diff = c - c_approx;
disp('The Difference between Coefficients:'); disp(coeff_diff);

%--- 2.e) 
% Compute the correlation matrix of the estimated coefficients

% Define the number of data points for approximation case
M = length(xi_approx);

% Calculate the estimated measurement error using the difference of A_approx*c_approx and y_approx
e = A_approx * c_approx - yi_approx;

% Calculate the variance using the formula e'*e/(M - N)
sigma_2 = dot(e,e)/(M - N);

% Calculate the correlation matrix R using the formula sigma_2*inv(A_approx'*A_approx)
R = sigma_2 * inv(A_approx' * A_approx);

% For comparation of correlation matrices
% (R, the correlation matrix of Approximation is already calculated above)
R_interp = sigma_2 * inv(A' * A);

disp('Correlation Matrix for Interpolation:'); disp(R_interp);
disp('Correlation Matrix for Approximation:'); disp(R);
disp('Difference between Inter. and Aprox. Correlation Matrices:'); disp(R -R_interp);

%--- 2.f)
N2 = 2; % Reduced N value to 2

% Calculate coefficients - LSE with new reduced N value
[c_N2, A_N2] = coefficients(xi_approx, yi_approx, N2);

% Evaluate the polynomial at x_value (18) using the dot product of powers_x_N2 and c_N2
powers_x_N2 = x_value.^(0:N2);
fx_18_approx_value_N2 = dot(powers_x_N2, c_N2);

% Calculate the estimated measurement error using the difference of A_N2*c_N2 and yi_approx
e_N2 = A_N2 * c_N2 - yi_approx;

% Calculate the variance using the formula e'*e/(M-N)
sigma2_N2 = dot(e_N2,e_N2)/(M-N2);

% Calculate the correlation matrix R2 using the formula sigma_2*inv(A_approx'*A_approx)
R2 = sigma2_N2 * inv(A_N2' * A_N2);

%% Graphically Comparation for coefficients ( 2.c) )
% Generate x values for plotting
x_values_interpolation = linspace(min(xi_interp), max(xi_interp), 100);
x_values_approximation = linspace(min(xi_approx), max(xi_approx), 100);

% Evaluate the polynomials for the generated x values
y_values_interpolation = polyval(c, x_values_interpolation);
y_values_approximation = polyval(c_approx, x_values_approximation);

% Plot the results
figure;
subplot(2, 1, 1);
plot(xi_interp, yi_interp, '*', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'DisplayName', 'Data Points (Interpolation)');
hold on;
plot(x_values_interpolation, y_values_interpolation, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Polynomial (Interpolation)');
xlabel('x'); ylabel('y'); title('Comparison of Polynomials (Interpolation)');
legend('show'); grid on; hold off;

subplot(2, 1, 2);
plot(xi_approx, yi_approx, '*', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'DisplayName', 'Data Points (Approximation)');
hold on;
plot(x_values_approximation, y_values_approximation, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Polynomial (Least-Squares Estimation)');
xlabel('x'); ylabel('y'); title('Comparison of Polynomials (Least-Squares Estimation)');
legend('show'); grid on; hold off;
sgtitle('Comparison of Lagrange Interpolation and Least-Squares Estimation Polynomials'); % Adjust layout


%% Which of both results are more reliable? (-> graphs, results etc. for 2.d) Interpolation vs Approximation)


%% Graphically comparing correlation matrices R (N=4) and R2 (N=2) for 2.f)
% Generate x values for plotting
generated_x_val = linspace(min(xi_approx), max(xi_approx), 100);

% Evaluate the polynomials for the generated x values
y_values_approx = polyval(c_approx, generated_x_val);
y_values_approx_reduced = polyval(c_N2, generated_x_val);

% Plot the results
figure;
plot(xi_approx, yi_approx, '*', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'DisplayName', 'Data Points');
hold on;
plot(generated_x_val, y_values_approx, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Polynomial (N = 4)');
plot(generated_x_val, y_values_approx_reduced, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Polynomial (N = 2)');
xlabel('x'); ylabel('y'); legend('show'); grid on; hold off;
sgtitle('Comparison of N=4 and N=2 in Approximation Case');

%%
% Graphical comparison of correlation matrices
figure;

% Plot correlation matrix for N=4
subplot(1,2,1);
imagesc(R);
title('Correlation Matrix for N=4');
colorbar;

% Plot correlation matrix for N=2 (reduced degree)
subplot(1,2,2);
imagesc(R2);
title('Correlation Matrix for N=2 (Reduced Degree)');
colorbar;

% Adjust figure layout
sgtitle('Graphical Comparison of Correlation Matrices');

% Optional: Display numerical values on each cell
disp('Correlation Matrix for N=4:');
disp(R);
disp('Correlation Matrix for N=2 (Reduced Degree):');
disp(R2);
