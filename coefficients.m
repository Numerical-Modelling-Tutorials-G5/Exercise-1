% This function calculates the coefficients of a polynomial that best fits
% the data points given by x and y, using the method of least squares.
% The function also returns the matrix A that represents the linear system
% of equations that relates the coefficients and the data points.
% The input arguments are:
% x - a vector of the x-coordinates of the data points
% y - a vector of the y-coordinates of the data points
% N - the maximum degree of the polynomial
% The output arguments are:
% c - a vector of the coefficients of the polynomial, in ascending order of power
% A - the matrix of the linear system of equations

function [c, A] = coefficients(x, y, N)
  
    % Creating the matrix A:
    A = vander(x); % Using the vander function to generate a Vandermonde matrix from x
    A = fliplr(A); % Flipping the matrix horizontally, so that the columns are in ascending order of power

    % Keeping only the first N+1 columns of A, corresponding to the powers from 0 to N
    A = A(:,1:N+1);
    % Computing the normal equation matrix A'A and its inverse
    A_norm = A' * A;
    A_norm_inv = inv(A_norm);  
    
    % Solving for the coefficients c using the formula c = (A'A)^-1 A'y
    c = A_norm_inv * A' * y;
end