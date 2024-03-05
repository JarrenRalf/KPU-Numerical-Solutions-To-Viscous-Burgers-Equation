function [L, U] = LUdecompCrout(A)
% This function LUdecompCrout, decomposes a given matrix, A, into a lower triangular matrix, L, and an
% upper triangular matrix, U. This function was modified and vectorized by Jarren Ralf.
%
%            A = Any N x M matrix with real number entries
%
% @author https://en.wikipedia.org/wiki/Crout_matrix_decomposition
        
[R, C] = size(A); % Retrieve the number of rows and columns of the given matrix
L = zeros(R, C);  % Initialize the matrix L with the same dimensions as matrix A
U = eye(R, C);    % Initialize the matrix U with the identity matrix

L(:, 1) = A(:, 1); % Extract the first column of A and set it as the first column of L
U(1, 2:end) = A(1, 2:end)/L(1, 1); % Set the fist row of the U matrix starting from the second column

for i = 2:R
    
    % Construct the lower triangular matrix, L
    for j = 2:i
        L(i, j) = A(i, j) - L(i, 1:j - 1)*U(1:j - 1, j);
    end

    % Construct the upper triangular matrix, U
    for j = i + 1:R
        U(i, j) = (A(i, j) - L(i, 1:i - 1)*U(1:i - 1, j))/L(i, i);
    end
end