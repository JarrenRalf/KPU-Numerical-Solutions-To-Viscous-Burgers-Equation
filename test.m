%%%%%%%%%%%%%%%%%%%%%%
%% Analytical Solution
%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

a  = 1;    % Amplitude of the sin initial condition
L  = 1;    % Length of the x-interval -- [0, L]
T  = 1;    % Length of the t-intercal -- [0, T]
c  = 1;    % Diffusion/Viscosity Constant
dt = 0.02; % This is delta t -- The size of the sub-interval in time
dx = 0.01; % This is delta x -- The size of the sub-interval in space
numTerms  = 100; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition
initialCond = @(x) a*sin(pi*x/L);
%initialCond = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t] = BurgersEq1D(f, c, L, T, dx, dt, numTerms); 

% Plot the surface
surf(t, x', u);
title(['The Solution of Burgers'' equation for \nu = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analytical Solution - Using Numerical Integration while performing the Cole-Hopf transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

a  = 1;    % Amplitude of the sin initial condition
L  = 1;    % Length of the x-interval -- [0, L]
T  = 1;    % Length of the t-intercal -- [0, T]
c  = 1;    % Diffusion/Viscosity Constant
dt = 0.02; % This is delta t -- The size of the sub-interval in time
dx = 0.01; % This is delta x -- The size of the sub-interval in space

% Number of terms in the finite Fourier series -- Errors occur when too many terms are used
numTerms  = 100; 

% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) a*sin(pi*z/L);
%u0 = @(z) 4*z.*(1 - z);

% Apply the Cole-Hopf Transformation to the initial condition
x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

% Solve the equation
[u, x, t] = BurgersEq1D_Numerical(x0, c, L, T, dx, dt, numTerms);

% Plot the surface
surf(t, x', u);
title(['The Solution of Burgers'' equation for \nu = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Animation of the Analytical Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

a  = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 1.5;   % Length of the t-intercal -- [0, T]
c  = 1;     % Diffusion/Viscosity Constant
dt = 0.002; % This is delta t -- The size of the sub-interval in time
dx = 0.001; % This is delta x -- The size of the sub-interval in space
numTerms  = 100; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) a*sin(pi*z/L);
%u0 = @(z) 4*z.*(1 - z);

% Apply the Cole-Hopf Transformation to the initial condition
x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

% Solve the equation
[u, x, t] = BurgersEq1D_Numerical(x0, c, L, T, dx, dt, numTerms);  

% Created an animated solution and save it as .gif
AnimatedPDESoln(u, u0, c, x, t, L, a, 'BurgersEquation.gif');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finite Elements Scheme -- Crout method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

a  = 1;      % Amplitude of the sin initial condition
L  = 1;      % Length of the x-interval -- [0, L]
T  = 1;      % Length of the t-intercal -- [0, T]
c  = 1;      % Diffusion/Viscosity Constant
dt = 0.0016; % This is delta t -- The size of the sub-interval in time
dx = 0.001;  % This is delta x -- The size of the sub-interval in space

% Chosoe the method     0 -> Explicit;     .5 -> Crank-Nicolson;     1 -> Implicit
theta = 0;

% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) a*sin(pi*z/L);
%u0 = @(z) 4*z.*(1 - z);

% Apply the Cole-Hopf Transformation to the initial condition
x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

% Solve the Heat equation via a finite elements scheme using the Crout method
[v, x, t] = HeatEq1D_FiniteElements_Crout(x0, c, L, T, dx, dt, theta);

% Compute the solution to Burgers' equation by applying the Cole-Hopf transformation
u = ColeHopfTransformation_Discrete(v, c, L, T, dx, dt);

% Plot the surface
surf(t, x', u);
title(['The Solution of Burgers'' equation for \nu = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')




%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finite Elements Scheme
%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

a  = 1;      % Amplitude of the sin initial condition
L  = 1;      % Length of the x-interval -- [0, L]
T  = 1;      % Length of the t-intercal -- [0, T]
c  = 1;      % Diffusion/Viscosity Constant
dt = 0.0016; % This is delta t -- The size of the sub-interval in time
dx = 0.01;   % This is delta x -- The size of the sub-interval in space

% Chosoe the method     0 -> Explicit;     .5 -> Crank-Nicolson;     1 -> Implicit
theta = 1;

% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) a*sin(pi*z/L);
%u0 = @(z) 4*z.*(1 - z);

% Apply the Cole-Hopf Transformation to the initial condition
x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

% Solve the Heat equation via a finite elements scheme
[v, x, t] = HeatEq1D_FiniteElements(x0, c, L, T, dx, dt, theta);

% Compute the solution to Burgers' equation by applying the Cole-Hopf transformation
u = ColeHopfTransformation_Discrete(v, c, L, T, dx, dt);

% Plot the surface
surf(t, x', u);
title(['The Solution of Burgers'' equation for \nu = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Crank-Nicolson Finite Difference Scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

a  = 1;    % Amplitude of the sin initial condition
L  = 1;    % Length of the x-interval -- [0, L]
T  = 2;    % Length of the t-intercal -- [0, T]
c  = 1/5;  % Diffusion/Viscosity Constant
dt = 0.05; % This is delta t -- The size of the sub-interval in time
dx = 0.02; % This is delta x -- The size of the sub-interval in space

% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) a*sin(pi*z/L);
%u0 = @(z) 4*z.*(1 - z);

% Apply the Cole-Hopf Transformation to the initial condition
x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

% Solve the heat equation using finite differences via the Crank-Nicolson scheme
[v, x, t] = HeatEq1D_CrankNicolson(x0, c, L, T, dx, dt);

% Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
u = ColeHopfTransformation_Discrete(v, c, L, T, dx, dt);

% Plot the surface
surf(t, x', u);
title(['The Solution of Burgers'' equation for \nu = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')