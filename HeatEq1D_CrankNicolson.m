function [v, x, t, N, M] = HeatEq1D_CrankNicolson(x0, c, L, T, dx, dt)
% This function solveHeatEq1D_CrankNicolson, solves the one dimensional Heat equation using the finite
% difference technique via the Crank-Nicolson scheme. This function outputs the solution to the heat
% equation as a matrix of points. In addition the user can also receive the x and t-values that the 
% solution is evaluated at. Lastly, the number of sub-intervals in the x and t dimensions can also be
% passed to the user.
%
%         x0 = The initial condition passed as vector -- positive real numbers
%         c  = The diffusion/viscosity constant       -- positive real number
%         L  = The end point of the interval in space -- positive real number -- [0, L]
%         T  = The end point of the interval in time  -- positive real number -- [0, T]
%         dx = The small change in x on the grid -- Delta x -- positive real number *small i.e. < L
%         dt = The small change in t on the grid -- Delta t -- positive real number *small i.e. < T
%
% @author Jarren Ralf

t = 0:dt:T;      % t-grid points
x = 0:dx:L;      % x-grid points
M = floor(T/dt); % The number of subintervals in time
N = floor(L/dx); % The number of subintervals in space
r = dt*c/dx^2;   % Defined constant used in the Crank-Nicolson scheme

% - MATLAB indices start at 1, NOT 0. Hence the code will have many '+1's for the purpose of re-indexing -
% Initialize the solution matrix
v = zeros(N + 1, M + 1);

% Build the tridiagonal matrix
mtx = full(gallery('tridiag', N + 1, -1/2*r, 1 + r, -1/2*r));
mtx(1    , 2) = -1*r;
mtx(N + 1, N) = -1*r;

% Use the initial condition to populate the first column of the solution matrix for t = 0
v(:, 0 + 1) = x0;

% Initialize the vector b
b = zeros(N + 1, 1);

% Apply the Crank-Nicolson Scheme
for j = 0:M - 1
    
    % The i = 0 case
    b(0 + 1) = r*v((0 + 1) + 1, j + 1) + (1 - r)*v(0 + 1, j + 1);

    % The i = 1, 2, ..., N - 1 case
    b(2:end-1) = r/2*v((2:end-1) + 1, j + 1) + (1 - r)*v(2:end-1, j + 1) + r/2*v((2:end-1) - 1, j + 1);
    
    % The i = N case
    b(N + 1) = (1 - r)*v(N + 1, j + 1) + r*v((N + 1) - 1, j + 1);

    % Solve the system for the next time step
    v(:, (j + 1) + 1) = mtx\b;
end