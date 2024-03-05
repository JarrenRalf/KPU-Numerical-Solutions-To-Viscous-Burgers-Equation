function x0 = ColeHopfTransformation_Numerical(u0, c, L, dx)
% This function coleHopfTransformation_Numerical, applies the Cole-Hopf transformation to a given 
% anonymous function for a chosen value of the dissipation coefficient, interval length and space step. 
% This function receives the initial condition for Burgers' equation and converts it into the 
% corresponding initial condition for the Heat equation. This method uses ode45 which is a high order 
% numerical ode solver.
%
%               u0 = The initial condition passed as an anonymous function
%               c  = The diffusion/viscosity constant       -- positive real number
%               L  = The end point of the interval in space -- positive real number -- [0, L]
%               dx = The small change in x on the grid -- Delta x -- positive real number *small i.e. < L
%
% @author Jarren Ralf

x = 0:dx:L; % x-grid points
y0 = 1;     % The iniital condition used to solve the ODE

% Rearrange the ODE in the form y' = f(z, y(z)) and set up the right hand side as an anonymous function
f = @(z, y) -(y*u0(z))/(2*c);

% Solve the ODE and output the solution as a column vector of of real numbers
[z, x0] = ode45(f, x, y0);