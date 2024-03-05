function u = ColeHopfTransformation_Discrete(v, c, L, T, dx, dt)
% This function coleHopfTransformation_Discrete, applies the discrete Cole-Hopf transformation to the 
% solution of the Diffusion equation, v. The output is the approximate solution to Burgers' equation.
%
%         v  = The solution to the Diffusion equation -- matrix of real numbers
%         c  = The diffusion/viscosity constant       -- positive real number
%         L  = The end point of the interval in space -- positive real number -- [0, L]
%         T  = The end point of the interval in time  -- positive real number -- [0, T]
%         dx = The small change in x on the grid -- Delta x -- positive real number *small i.e. < L
%         dt = The small change in t on the grid -- Delta t -- positive real number *small i.e. < T
%
% @author Jarren Ralf

M  = floor(T/dt);  % The number of subintervals in time
N  = floor(L/dx);  % The number of subintervals in space

u = zeros(N + 1, M + 1); % Initialize the solution matrix for Burgers' equation

% Apply the cole-hopf transformation via a forward difference approximation for the first derivative
u(0 + 1, :) = -2*c*((v(1 + 1, :) - v(0 + 1, :))./(dx*v(0 + 1, :)));

% Apply the cole-hopf transformation via a centered difference approximation for the first derivative
for i = 1:N - 1
    u(i + 1, :) = -1*c*((v((i + 1) + 1, :) - v((i + 1) - 1, :))./(dx*v(i + 1, :)));
end

% Apply the cole-hopf transformation via a backward difference approximation for the first derivative
u(N + 1, :) = -2*c*((v(N + 1, :) - v((N - 1) + 1, :))./(dx*v(N + 1, :)));