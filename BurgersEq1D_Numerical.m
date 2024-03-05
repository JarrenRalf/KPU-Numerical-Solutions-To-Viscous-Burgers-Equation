function [u, x, t, denominator] = BurgersEq1D_Numerical(x0, c, L, T, dx, dt, N)
% This function solnBurgersEq1D, solves the one-dimentional viscous Burgers' equation using a finite 
% Fourier series approximatation to the diffusion equation. This function outputs the solution to the heat
% equation as a matrix of points. In addition the user can also receive the x and t-values that the 
% solution is evaluated at. Lastly, the user could be passed the solution to the diffusion equation,
% denoted here by 'denominator'.
%
%         x0 = The initial condition passed as vector -- positive real numbers
%         c  = The diffusion/viscosity constant       -- positive real number
%         L  = The end point of the interval in space -- positive real number -- [0, L]
%         T  = The end point of the interval in time  -- positive real number -- [0, T]
%         dx = The small change in x on the grid -- Delta x -- positive real number *small i.e. < L
%         dt = The small change in t on the grid -- Delta t -- positive real number *small i.e. < T
%         N  = The number of terms used in the Fourier series approximation -- large positive real number
%
% @author Jarren Ralf

t = 0:dt:T; % t-grid points
x = 0:dx:L; % x-grid points
x = x';         % transpose to make column vector

a0 = trapz(x, x0)/L; % Calculate the a0 constant

% Initialize the coefficient, numerator, and denominator matrices
          a =   zeros(1, N);
  numerator =   zeros(length(x), length(t));
denominator = a0*ones(length(x), length(t)); 

for n = 1:N

    % Calculate the appropriate Fourier coefficient at step n
    a(n) = 2*trapz(x, x0.*cos(n*pi*x/L))/L;

    % Update the grid points, each node in the mesh, of the numerator/denominator 
    % matrix with the next function evaluation of the Fourier series at step n 
      numerator =   numerator + n*a(n)*exp(-1*n^2*pi^2*c*t/L^2).*sin(n*pi*x/L);
    denominator = denominator +   a(n)*exp(-1*n^2*pi^2*c*t/L^2).*cos(n*pi*x/L);
end

% Construct the approximation of the analytical solution of Burgers' equation
u = 2*pi*c*numerator/L./denominator;