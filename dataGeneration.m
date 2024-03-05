%% Table 1
clear; clc; close all;

a  = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 0.1;   % Length of the t-intercal -- [0, T]
c  = 1;     % Diffusion/Viscosity Constant
dx = 1/10;  % This is delta x -- The size of the sub-interval in space
dt = 0.001; % This is delta t -- The size of the sub-interval in time
Col = 5;    % Set the number of columns in the table
ROW = 9;    % This is the number of rows -- this value is treated like a constant * Do not change
numTerms  = 100; % Number of terms in the finite Fourier series

% Initialize the table
table = zeros(ROW, Col);
    
% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) a*sin(pi*z/L);

for j = 1:Col
    
    % Apply the Cole-Hopf Transformation to the initial condition
    x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

    % Solve the heat equation using finite differences via the Crank-Nicolson scheme
    [v, x, t, N] = HeatEq1D_CrankNicolson(x0, c, L, T, dx, dt); 

    % Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
    u = BurgersEq1D_Numerical(x0, c, L, T, dx, dt, numTerms);
    u_approx = ColeHopfTransformation_Discrete(v, c, L, T, dx, dt);

    % Calculate the approprate abscissae values
    abscissae = N/(ROW + 1) + 1:N/(ROW + 1):N;
    
    % Build the table
    table(:, j) = abs(u_approx(abscissae', end) - u(abscissae', end));
    
    % Update delta x
    dx = dx/2;
end

% Display the table
disp(num2str(table, '%15.6f'));
%% Table 2
clear; clc; close all;

a  = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 0.1;   % Length of the t-intercal -- [0, T]
c  = 1/15;  % Diffusion/Viscosity Constant
dx = 1/10;  % This is delta x -- The size of the sub-interval in space
dt = 0.001; % This is delta t -- The size of the sub-interval in time
Col = 5;    % Set the number of columns in the table
ROW = 9;    % This is the number of rows -- this value is treated like a constant * Do not change
numTerms  = 100; % Number of terms in the finite Fourier series

% Initialize the table
table = zeros(ROW, Col);
    
% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) a*sin(pi*z/L);

for j = 1:Col
    
    % Apply the Cole-Hopf Transformation to the initial condition
    x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

    % Solve the heat equation using finite differences via the Crank-Nicolson scheme
    [v, x, t, N] = HeatEq1D_CrankNicolson(x0, c, L, T, dx, dt); 

    % Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
    u = BurgersEq1D_Numerical(x0, c, L, T, dx, dt, numTerms);
    u_approx = ColeHopfTransformation_Discrete(v, c, L, T, dx, dt);

    % Calculate the approprate abscissae values
    abscissae = N/(ROW + 1) + 1:N/(ROW + 1):N;
    
    % Build the table
    table(:, j) = abs(u_approx(abscissae', end) - u(abscissae', end));
    
    % Update delta x
    dx = dx/2;
end

% Display the table
disp(num2str(table, '%15.6f'));
%% Table 3
clear; clc; close all;

a  = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 0.01;  % Length of the t-intercal -- [0, T]
c  = 1/15;  % Diffusion/Viscosity Constant
dx = 1/10;  % This is delta x -- The size of the sub-interval in space
dt = 0.001; % This is delta t -- The size of the sub-interval in time
Col = 5;    % Set the number of columns in the table
ROW = 9;    % This is the number of rows -- this value is treated like a constant * Do not change
numTerms  = 100; % Number of terms in the finite Fourier series

% Initialize the table
table = zeros(ROW, Col);
    
% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) a*sin(pi*z/L);

for j = 1:Col
    
    % Apply the Cole-Hopf Transformation to the initial condition
    x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

    % Solve the heat equation using finite differences via the Crank-Nicolson scheme
    [v, x, t, N] = HeatEq1D_CrankNicolson(x0, c, L, T, dx, dt); 

    % Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
    u = BurgersEq1D_Numerical(x0, c, L, T, dx, dt, numTerms);
    u_approx = ColeHopfTransformation_Discrete(v, c, L, T, dx, dt);

    % Calculate the approprate abscissae values
    abscissae = N/(ROW + 1) + 1:N/(ROW + 1):N;
    
    % Build the table
    table(:, j) = abs(u_approx(abscissae', end) - u(abscissae', end));
    
    % Update delta x
    dx = dx/2;
end

% Display the table
disp(num2str(table, '%15.6f'));
%% Table 4
clear; clc; close all;

a  = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 0.1;   % Length of the t-intercal -- [0, T]
c  = 1;     % Diffusion/Viscosity Constant
dx = 1/10;  % This is delta x -- The size of the sub-interval in space
dt = 0.001; % This is delta t -- The size of the sub-interval in time
Col = 5;    % Set the number of columns in the table
ROW = 9;    % This is the number of rows -- this value is treated like a constant * Do not change
numTerms  = 100; % Number of terms in the finite Fourier series

% Initialize the table
table = zeros(ROW, Col);
    
% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) 4*z.*(1 - z);

for j = 1:Col
    
    % Apply the Cole-Hopf Transformation to the initial condition
    x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

    % Solve the heat equation using finite differences via the Crank-Nicolson scheme
    [v, x, t, N] = HeatEq1D_CrankNicolson(x0, c, L, T, dx, dt); 

    % Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
    u = BurgersEq1D_Numerical(x0, c, L, T, dx, dt, numTerms);
    u_approx = ColeHopfTransformation_Discrete(v, c, L, T, dx, dt);

    % Calculate the approprate abscissae values
    abscissae = N/(ROW + 1) + 1:N/(ROW + 1):N;
    
    % Build the table
    table(:, j) = abs(u_approx(abscissae', end) - u(abscissae', end));
    
    % Update delta x
    dx = dx/2;
end

% Display the table
disp(num2str(table, '%15.6f'));
%% Table 5
clear; clc; close all;

a  = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 0.1;   % Length of the t-intercal -- [0, T]
c  = 1/15;  % Diffusion/Viscosity Constant
dx = 1/10;  % This is delta x -- The size of the sub-interval in space
dt = 0.001; % This is delta t -- The size of the sub-interval in time
Col = 5;    % Set the number of columns in the table
ROW = 9;    % This is the number of rows -- this value is treated like a constant * Do not change
numTerms  = 100; % Number of terms in the finite Fourier series

% Initialize the table
table = zeros(ROW, Col);
    
% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) 4*z.*(1 - z);

for j = 1:Col
    
    % Apply the Cole-Hopf Transformation to the initial condition
    x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

    % Solve the heat equation using finite differences via the Crank-Nicolson scheme
    [v, x, t, N] = HeatEq1D_CrankNicolson(x0, c, L, T, dx, dt); 

    % Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
    u = BurgersEq1D_Numerical(x0, c, L, T, dx, dt, numTerms);
    u_approx = ColeHopfTransformation_Discrete(v, c, L, T, dx, dt);

    % Calculate the approprate abscissae values
    abscissae = N/(ROW + 1) + 1:N/(ROW + 1):N;
    
    % Build the table
    table(:, j) = abs(u_approx(abscissae', end) - u(abscissae', end));
    
    % Update delta x
    dx = dx/2;
end

% Display the table
disp(num2str(table, '%15.6f'));
%% Table 6
clear; clc; close all;

a  = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 0.01;  % Length of the t-intercal -- [0, T]
c  = 1/15;  % Diffusion/Viscosity Constant
dx = 1/10;  % This is delta x -- The size of the sub-interval in space
dt = 0.001; % This is delta t -- The size of the sub-interval in time
Col = 5;    % Set the number of columns in the table
ROW = 9;    % This is the number of rows -- this value is treated like a constant * Do not change
numTerms  = 100; % Number of terms in the finite Fourier series

% Initialize the table
table = zeros(ROW, Col);
    
% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) 4*z.*(1 - z);

for j = 1:Col
    
    % Apply the Cole-Hopf Transformation to the initial condition
    x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

    % Solve the heat equation using finite differences via the Crank-Nicolson scheme
    [v, x, t, N] = HeatEq1D_CrankNicolson(x0, c, L, T, dx, dt); 

    % Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
    u = BurgersEq1D_Numerical(x0, c, L, T, dx, dt, numTerms);
    u_approx = ColeHopfTransformation_Discrete(v, c, L, T, dx, dt);

    % Calculate the approprate abscissae values
    abscissae = N/(ROW + 1) + 1:N/(ROW + 1):N;
    
    % Build the table
    table(:, j) = abs(u_approx(abscissae', end) - u(abscissae', end));
    
    % Update delta x
    dx = dx/2;
end

% Display the table
disp(num2str(table, '%15.6f'));
%% Table 7
clear; clc; close all;

a  = 1;      % Amplitude of the sin initial condition
L  = 1;      % Length of the x-interval -- [0, L]
T  = 3;      % Length of the t-intercal -- [0, T]
c  = 1;      % Diffusion/Viscosity Constant
dx = 0.0125; % This is delta x -- The size of the sub-interval in space
dt = 0.01;   % This is delta t -- The size of the sub-interval in time
Col = 5;     % Set the number of columns in the table
ROW = 9;     % This is the number of rows -- this value is treated like a constant * Do not change
numTerms  = 100; % Number of terms in the finite Fourier series

% Initialize the table
table = zeros(ROW, Col);
    
% Toggle the chosen initial condition (feel free to define your own!)
u0 = @(z) a*sin(pi*z/L);

for j = 1:Col
    
    % Apply the Cole-Hopf Transformation to the initial condition
    x0 = ColeHopfTransformation_Numerical(u0, c, L, dx);

    % Solve the heat equation using finite differences via the Crank-Nicolson scheme
    [v, x, t, N] = HeatEq1D_CrankNicolson(x0, c, L, T, dx, dt); 

    % Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
    u = BurgersEq1D_Numerical(x0, c, L, T, dx, dt, numTerms);
    u_approx = ColeHopfTransformation_Discrete(v, c, L, T, dx, dt);

    % Calculate the approprate abscissae values
    abscissae = N/(ROW + 1) + 1:N/(ROW + 1):N;
    
    % Build the table
    table(:, j) = abs(u_approx(abscissae', end) - u(abscissae', end));
    
    % Update delta x
    dx = dx/2;
end

% Display the table
disp(num2str(table, '%15.6f'));