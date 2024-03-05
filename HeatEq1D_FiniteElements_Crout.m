function [v, x, t, N, M] = HeatEq1D_FiniteElements_Crout(x0, c, L, T, dx, dt, theta)
% This function solveHeatEq1D_FiniteElements_Crout, solves the one dimensional Heat equation using 
% a finite elements technique which uses the Crout method for solving the system of linear equations. 
% This function outputs the solution to the heat equation as a matrix of points. In addition the user 
% can also receive the x and t-values that the solution is evaluated at. Lastly, the number of 
% sub-intervals in the x and t dimensions can also be passed to the user.
%
%         x0 = The initial condition passed as vector -- positive real numbers
%         c  = The diffusion/viscosity constant       -- positive real number
%         L  = The end point of the interval in space -- positive real number -- [0, L]
%         T  = The end point of the interval in time  -- positive real number -- [0, T]
%         dx = The small change in x on the grid -- Delta x -- positive real number *small i.e. < L
%         dt = The small change in t on the grid -- Delta t -- positive real number *small i.e. < T
%      theta = The choice of weight factor -- positive real number in [0, 1]
%
% @author Jarren Ralf

if theta < 0.5
    if (c*dt)/dx^2 > (1 - 2*theta)/6
        disp('The finite element equation is unstable!');
    end
end

x = 0:dx:L;      % x-grid points
t = 0:dt:T;      % x-grid points
N = floor(L/dx); % The number of subintervals in space
M = floor(T/dt); % The number of subintervals in time

% Construct the matrices C and K (see my report for derivation of these quantities)
cii =  2*dx/3;         
cij =    dx/6;
kii =  2/dx;
kij = -1/dx;
  C = full(gallery('tridiag', N + 1, cij, cii, cij));
  K = full(gallery('tridiag', N + 1, kij, kii, kij));

% Construct the matrix D and E matrix
D = C + c*     theta *dt*K;
E = C - c*(1 - theta)*dt*K;

[L, U] = LUdecompCrout(D); % Build the Lower and Upper tridiagonal matrices

v = zeros(N + 1, M + 1); % Initialize the solution matrix
y = zeros(N + 1,     1); % Initialize the y vector

v(:, 1) = x0; % Use the initial condition to populate the first column of the solution matrix for t = 0

% Loop through each time step
for j = 1:M
      b  = E*v(:, j);    % Initialize the b matrix
    y(1) = b(1)/L(1, 1); % Set the first y-value i.e. for n = 1
    
    % Loop through the rest of the space steps to build the y vector
    for i = 2:N + 1
        y(i) = (b(i) - L(i, i - 1)*y(i - 1))/L(i, i);
    end
    
    v(N + 1, j + 1) = y(N + 1); % Set the last x-value i.e. for n = N + 1
    
    % Loop backwards through the rest of the space steps to build the x vector
    for i = N:-1:1
        v(i, j + 1) = y(i) - U(i, i + 1)*v(i + 1, j + 1);
    end
end