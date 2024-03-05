%% First poster image 
clear; clc; close all;

L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = .0543; % Diffusion/Viscosity Constant
dt = 0.25;  % This is delta t -- The size of the sub-interval in time
dx = 0.01;  % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);  

% Set initial condition
y = initialCond(x);

% Retreive the number of columns of the solution matrix
numColumns = size(u, 2); 

% Initialize the counter variable which chooses the subplots
j = 1;

for i = 1:numColumns
    
    % Choose only these four images
    if t(i) == 0.0 || t(i) == 0.25 || t(i) == 0.50 || t(i) == 2.0
        
        % Pass the updated function values to y, i.e. the solution changing in time
        y = u(:, i); 

        % Create a subplot
        subplot(1, 4, j);
        
        % Updates the data source, i.e. the y-values from above and label the axis
        plot(x, y);
        axis([0 L 0 1]);
        refreshdata
        sgtitle(['The Solution of Burgers'' equation for \color{red}\kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
        title(['t = ', num2str(t(i), '%5.3f')])
        xlabel('x')
        ylabel('u(x, t)')
    
        % Increment the counter of subplots
        j = j + 1;
    end
end
%% Second poster image
clear; clc; close all;

L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = .0544; % Diffusion/Viscosity Constant
dt = 0.25;  % This is delta t -- The size of the sub-interval in time
dx = 0.01;  % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t] = BurgersEq1D(f, c, L, T, dx, dt, numTerms); 

% Plot the surface
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
%% Image # 1 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1;     % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 1a for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/2;   % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 1b for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/3;   % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 1c for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/4;   % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 2 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/5;   % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 3 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/10;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 4 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/15;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 5 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/20;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 6 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/25;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 7 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/30;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 8 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/35;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 9 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/40;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 10 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/45;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 11 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/50;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 12 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/55;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 13 for presentation -- Different Initial Conditions
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/60;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond1 = @(x) u0*sin(pi*x/L);
initialCond2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
f1 = ColeHopfTransformation(initialCond1, c);
f2 = ColeHopfTransformation(initialCond2, c);

% Solve the equation
[u1, x1, t1] = BurgersEq1D(f1, c, L, T, dx, dt, numTerms);
[u2, x2, t2] = BurgersEq1D(f2, c, L, T, dx, dt, numTerms); 

% Plot the surface
subplot(1, 2, 1)
surf(t1, x1', u1);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t2, x2', u2);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(80, 10)
axis([0 T 0 L 0 u0])
%% Image # 14 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1;     % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 15 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/2;   % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 16 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/3;   % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 17 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/4;   % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 18 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/5;   % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 19 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/10;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 20 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/15;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 21 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/20;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 22 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/25;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 23 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1;     % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Solve the equation
[u2, x2, t2, denominator2] = BurgersEq1D(initialCond, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator2);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 24 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/5;   % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Solve the equation
[u2, x2, t2, denominator2] = BurgersEq1D(initialCond, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator2);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 25 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/10;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation *** set t = 0.001:dt:T; and x = 0.001:dx:L;
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Solve the equation
[u2, x2, t2, denominator2] = BurgersEq1D(initialCond, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator2);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 26 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/15;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation *** set t = 0.001:dt:T; and x = 0.001:dx:L;
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Solve the equation
[u2, x2, t2, denominator2] = BurgersEq1D(initialCond, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator2);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 27 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/20;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation *** set t = 0.001:dt:T; and x = 0.001:dx:L;
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Solve the equation
[u2, x2, t2, denominator2] = BurgersEq1D(initialCond, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator2);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 28 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/25;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation *** set t = 0.001:dt:T; and x = 0.001:dx:L;
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Solve the equation
[u2, x2, t2, denominator2] = BurgersEq1D(initialCond, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator2);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 29 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/30;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation *** set t = 0.001:dt:T; and x = 0.001:dx:L;
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Solve the equation
[u2, x2, t2, denominator2] = BurgersEq1D(initialCond, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator2);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 30 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/35;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation *** set t = 0.001:dt:T; and x = 0.001:dx:L;
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Solve the equation
[u2, x2, t2, denominator2] = BurgersEq1D(initialCond, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator2);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 31 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/40;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation *** set t = 0.001:dt:T; and x = 0.001:dx:L;
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Solve the equation
[u2, x2, t2, denominator2] = BurgersEq1D(initialCond, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator2);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 32 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/45;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation *** set t = 0.001:dt:T; and x = 0.001:dx:L;
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Solve the equation
[u2, x2, t2, denominator2] = BurgersEq1D(initialCond, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator2);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 33 for presentation -- Burgers' vs. Diffusion
clear; clc; close all;

u0 = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 4;     % Length of the t-intercal -- [0, T]
c  = 1/50;  % Diffusion/Viscosity Constant
dt = 0.05;  % This is delta t -- The size of the sub-interval in time
dx = 0.015; % This is delta x -- The size of the sub-interval in space
numTerms  = 1000; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
initialCond = @(x) u0*sin(pi*x/L);

% Apply the Cole-Hopf Transformation to the initial condition
f = ColeHopfTransformation(initialCond, c);

% Solve the equation *** set t = 0.001:dt:T; and x = 0.001:dx:L;
[u, x, t, denominator] = BurgersEq1D(f, c, L, T, dx, dt, numTerms);

% Solve the equation
[u2, x2, t2, denominator2] = BurgersEq1D(initialCond, c, L, T, dx, dt, numTerms);

% Plot the surface
subplot(1, 2, 1)
surf(t, x', u);
title(['The Solution of Burgers'' equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('u(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])

subplot(1, 2, 2)
surf(t, x', denominator2);
title(['The Solution of the Diffusion equation for \kappa = ', num2str(c, '%f'), ' (\itRe = ', num2str(1/c), ')'])
xlabel('t')
ylabel('x')
zlabel('v(x, t)')
view(130, 10)
axis([0 T 0 L 0 u0])
%% Image # 1 for thesis -- Burgers' vs. Diffusion
clear; clc; close all;

a  = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 1;     % Length of the t-intercal -- [0, T]
c  = 1;     % Diffusion/Viscosity Constant
dt = 0.001; % This is delta t -- The size of the sub-interval in time
dx = 0.001; % This is delta x -- The size of the sub-interval in space
numTerms  = 100; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
u0_1 = @(x) a*sin(pi*x/L);
u0_2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
x1 = ColeHopfTransformation_Numerical(u0_1, c, L, dx);
x2 = ColeHopfTransformation_Numerical(u0_2, c, L, dx);

% Solve the equation
 u1        = BurgersEq1D_Numerical(x1, c, L, T, dx, dt, numTerms);  
[u2, x, t] = BurgersEq1D_Numerical(x2, c, L, T, dx, dt, numTerms);

% Retreive the number of columns of the solution matrix
numColumns = size(u1, 2); 

% Initialize the counter variable which chooses the subplots
j = 1;

for i = 1:numColumns

    % Choose only these four images
    if t(i) == 0.001 || t(i) == 0.05 || t(i) == 0.1 || t(i) == 0.40 

        % Pass the updated function values to y, i.e. the solution changing in time
        y1 = u1(:, i);
        y2 = u2(:, i); 

        % Updates the data source, i.e. the y-values from above and label the axis
        plot1 = subplot(1, 2, 1);
        txt = ['t = ', num2str(t(i), '%0.3f')];
        plot(x, y1, 'DisplayName',txt);
        axis([0 L 0 1]);
        title(['The Solution of Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
        xlabel('x')
        ylabel('u(x, t)')

        plot2 = subplot(1, 2, 2);
        txt = ['t = ', num2str(t(i), '%0.3f')];
        plot(x, y2, 'DisplayName',txt);
        axis([0 L 0 1]);
        title(['The Solution of Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
        xlabel('x')
        ylabel('u(x, t)')
        
        hold(plot1, 'on')
        hold(plot2, 'on')

        % Increment the counter of subplots
        j = j + 1;
    end
    if j == 5
        break;
    end
end
    
hold(plot1, 'off')
lgd1 = legend(plot1, 'show');
lgd1.Location = 'northwest';

hold(plot2, 'off')
lgd2 = legend(plot2, 'show');
lgd2.Location = 'northwest';
%% Image # 2 for thesis -- Burgers' vs. Diffusion
clear; clc; close all;

a  = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 2;     % Length of the t-intercal -- [0, T]
c  = 1/25;     % Diffusion/Viscosity Constant
dt = 0.001; % This is delta t -- The size of the sub-interval in time
dx = 0.001; % This is delta x -- The size of the sub-interval in space
numTerms  = 100; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
u0_1 = @(x) a*sin(pi*x/L);
u0_2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
x1 = ColeHopfTransformation_Numerical(u0_1, c, L, dx);
x2 = ColeHopfTransformation_Numerical(u0_2, c, L, dx);

% Solve the equation
 u1        = BurgersEq1D_Numerical(x1, c, L, T, dx, dt, numTerms);  
[u2, x, t] = BurgersEq1D_Numerical(x2, c, L, T, dx, dt, numTerms);

% Retreive the number of columns of the solution matrix
numColumns = size(u1, 2); 

% Initialize the counter variable which chooses the subplots
j = 1;

for i = 1:numColumns

    % Choose only these four images
    if  t(i) == 0.05 || t(i) == 0.40 || t(i) == 1.0 || t(i) == 2.0

        % Pass the updated function values to y, i.e. the solution changing in time
        y1 = u1(:, i);
        y2 = u2(:, i); 

        % Updates the data source, i.e. the y-values from above and label the axis
        plot1 = subplot(1, 2, 1);
        txt = ['t = ', num2str(t(i), '%0.2f')];
        plot(x, y1, 'DisplayName',txt);
        axis([0 L 0 1]);
        title(['The Solution of Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
        xlabel('x')
        ylabel('u(x, t)')

        plot2 = subplot(1, 2, 2);
        txt = ['t = ', num2str(t(i), '%0.2f')];
        plot(x, y2, 'DisplayName',txt);
        axis([0 L 0 1]);
        title(['The Solution of Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
        xlabel('x')
        ylabel('u(x, t)')
        
        hold(plot1, 'on')
        hold(plot2, 'on')

        % Increment the counter of subplots
        j = j + 1;
    end
    if j == 5
        break;
    end
end
    
hold(plot1, 'off')
lgd1 = legend(plot1, 'show');
lgd1.Location = 'northwest';

hold(plot2, 'off')
lgd2 = legend(plot2, 'show');
lgd2.Location = 'northwest';
%% Image # 3 for thesis -- Burgers' vs. Diffusion
clear; clc; close all;

a  = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 2;     % Length of the t-intercal -- [0, T]
c  = 1/50;     % Diffusion/Viscosity Constant
dt = 0.001; % This is delta t -- The size of the sub-interval in time
dx = 0.001; % This is delta x -- The size of the sub-interval in space
numTerms  = 100; % Number of terms in the finite Fourier series

% Toggle the chosen initial condition (feel free to define your own!)
u0_1 = @(x) a*sin(pi*x/L);
u0_2 = @(x) 4*x.*(1 - x);

% Apply the Cole-Hopf Transformation to the initial condition
x1 = ColeHopfTransformation_Numerical(u0_1, c, L, dx);
x2 = ColeHopfTransformation_Numerical(u0_2, c, L, dx);

% Solve the equation
 u1        = BurgersEq1D_Numerical(x1, c, L, T, dx, dt, numTerms);  
[u2, x, t] = BurgersEq1D_Numerical(x2, c, L, T, dx, dt, numTerms);

% Retreive the number of columns of the solution matrix
numColumns = size(u1, 2); 

% Initialize the counter variable which chooses the subplots
j = 1;

for i = 1:numColumns

    % Choose only these four images
    if  t(i) == 0.05 || t(i) == 0.40 || t(i) == 1.0 || t(i) == 2.0

        % Pass the updated function values to y, i.e. the solution changing in time
        y1 = u1(:, i);
        y2 = u2(:, i); 

        % Updates the data source, i.e. the y-values from above and label the axis
        plot1 = subplot(1, 2, 1);
        txt = ['t = ', num2str(t(i), '%0.2f')];
        plot(x, y1, 'DisplayName',txt);
        axis([0 L 0 1]);
        title(['The Solution of Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
        xlabel('x')
        ylabel('u(x, t)')

        plot2 = subplot(1, 2, 2);
        txt = ['t = ', num2str(t(i), '%0.2f')];
        plot(x, y2, 'DisplayName',txt);
        axis([0 L 0 1]);
        title(['The Solution of Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
        xlabel('x')
        ylabel('u(x, t)')
        
        hold(plot1, 'on')
        hold(plot2, 'on')

        % Increment the counter of subplots
        j = j + 1;
    end
    if j == 5
        break;
    end
end
    
hold(plot1, 'off')
lgd1 = legend(plot1, 'show');
lgd1.Location = 'northwest';

hold(plot2, 'off')
lgd2 = legend(plot2, 'show');
lgd2.Location = 'northwest';
%% Image # 4 for thesis -- Error in Finite Difference Scheme
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
    
% Toggle the chosen initial condition (feel free to define your own!)
u0_1 = @(z) a*sin(pi*z/L);
u0_2 = @(z) 4*z.*(1 - z);

hold on

for j = 1:Col
    
    % Apply the Cole-Hopf Transformation to the initial condition
    x1 = ColeHopfTransformation_Numerical(u0_1, c, L, dx);
    x2 = ColeHopfTransformation_Numerical(u0_2, c, L, dx);

    % Solve the heat equation using finite differences via the Crank-Nicolson scheme
     v1           = HeatEq1D_CrankNicolson(x1, c, L, T, dx, dt);
    [v2, x, t, N] = HeatEq1D_CrankNicolson(x2, c, L, T, dx, dt);

    % Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
    u1 = BurgersEq1D_Numerical(x1, c, L, T, dx, dt, numTerms);
    u2 = BurgersEq1D_Numerical(x2, c, L, T, dx, dt, numTerms);
    u_approx1 = ColeHopfTransformation_Discrete(v1, c, L, T, dx, dt);
    u_approx2 = ColeHopfTransformation_Discrete(v2, c, L, T, dx, dt);

    % Calculate the approprate abscissae values
    abscissae = N/(ROW + 1) + 1:N/(ROW + 1):N;
    
    % Plot the error
    plot1 = subplot(1, 2, 1);
    txt = ['N = ', num2str(N)];
    plot(abs(u_approx1(abscissae', end) - u1(abscissae', end)), 'DisplayName', txt);
    title(['The Error in Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
    xlabel('x')
    xticks(1:9)
    xticklabels({'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'})
    ylabel('Absolute Error')
    
    plot2 = subplot(1, 2, 2);
    txt = ['N = ', num2str(N)];
    plot(abs(u_approx2(abscissae', end) - u2(abscissae', end)), 'DisplayName', txt);
    title(['The Error in Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
    xlabel('x')
    xticks(1:9)
    xticklabels({'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'})
    ylabel('Absolute Error')
    
    hold(plot1, 'on')
    hold(plot2, 'on')
    
    % Update delta x
    dx = dx/2;
end

hold(plot1, 'off')
lgd1 = legend(plot1, 'show');
lgd1.Location = 'northwest';

hold(plot2, 'off')
lgd2 = legend(plot2, 'show');
lgd2.Location = 'northwest';
%% Image # 5 for thesis -- Error in Finite Difference Scheme
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
    
% Toggle the chosen initial condition (feel free to define your own!)
u0_1 = @(z) a*sin(pi*z/L);
u0_2 = @(z) 4*z.*(1 - z);

hold on

for j = 1:Col
    
    % Apply the Cole-Hopf Transformation to the initial condition
    x1 = ColeHopfTransformation_Numerical(u0_1, c, L, dx);
    x2 = ColeHopfTransformation_Numerical(u0_2, c, L, dx);

    % Solve the heat equation using finite differences via the Crank-Nicolson scheme
     v1           = HeatEq1D_CrankNicolson(x1, c, L, T, dx, dt);
    [v2, x, t, N] = HeatEq1D_CrankNicolson(x2, c, L, T, dx, dt);

    % Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
    u1 = BurgersEq1D_Numerical(x1, c, L, T, dx, dt, numTerms);
    u2 = BurgersEq1D_Numerical(x2, c, L, T, dx, dt, numTerms);
    u_approx1 = ColeHopfTransformation_Discrete(v1, c, L, T, dx, dt);
    u_approx2 = ColeHopfTransformation_Discrete(v2, c, L, T, dx, dt);

    % Calculate the approprate abscissae values
    abscissae = N/(ROW + 1) + 1:N/(ROW + 1):N;
    
    % Plot the error
    plot1 = subplot(1, 2, 1);
    txt = ['N = ', num2str(N)];
    plot(abs(u_approx1(abscissae', end) - u1(abscissae', end)), 'DisplayName', txt);
    title(['The Error in Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
    xlabel('x')
    xticks(1:9)
    xticklabels({'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'})
    ylabel('Absolute Error')
    
    plot2 = subplot(1, 2, 2);
    txt = ['N = ', num2str(N)];
    plot(abs(u_approx2(abscissae', end) - u2(abscissae', end)), 'DisplayName', txt);
    title(['The Error in Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
    xlabel('x')
    xticks(1:9)
    xticklabels({'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'})
    ylabel('Absolute Error')
    
    hold(plot1, 'on')
    hold(plot2, 'on')
    
    % Update delta x
    dx = dx/2;
end

hold(plot1, 'off')
lgd1 = legend(plot1, 'show');
lgd1.Location = 'northwest';

hold(plot2, 'off')
lgd2 = legend(plot2, 'show');
lgd2.Location = 'northwest';
%% Image # 6 for thesis -- Error in Finite Difference Scheme
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
    
% Toggle the chosen initial condition (feel free to define your own!)
u0_1 = @(z) a*sin(pi*z/L);
u0_2 = @(z) 4*z.*(1 - z);

hold on

for j = 1:Col
    
    % Apply the Cole-Hopf Transformation to the initial condition
    x1 = ColeHopfTransformation_Numerical(u0_1, c, L, dx);
    x2 = ColeHopfTransformation_Numerical(u0_2, c, L, dx);

    % Solve the heat equation using finite differences via the Crank-Nicolson scheme
     v1           = HeatEq1D_CrankNicolson(x1, c, L, T, dx, dt);
    [v2, x, t, N] = HeatEq1D_CrankNicolson(x2, c, L, T, dx, dt);

    % Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
    u1 = BurgersEq1D_Numerical(x1, c, L, T, dx, dt, numTerms);
    u2 = BurgersEq1D_Numerical(x2, c, L, T, dx, dt, numTerms);
    u_approx1 = ColeHopfTransformation_Discrete(v1, c, L, T, dx, dt);
    u_approx2 = ColeHopfTransformation_Discrete(v2, c, L, T, dx, dt);

    % Calculate the approprate abscissae values
    abscissae = N/(ROW + 1) + 1:N/(ROW + 1):N;
    
    % Plot the error
    plot1 = subplot(1, 2, 1);
    txt = ['N = ', num2str(N)];
    plot(abs(u_approx1(abscissae', end) - u1(abscissae', end)), 'DisplayName', txt);
    title(['The Error in Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
    xlabel('x')
    xticks(1:9)
    xticklabels({'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'})
    ylabel('Absolute Error')
    
    plot2 = subplot(1, 2, 2);
    txt = ['N = ', num2str(N)];
    plot(abs(u_approx2(abscissae', end) - u2(abscissae', end)), 'DisplayName', txt);
    title(['The Error in Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
    xlabel('x')
    xticks(1:9)
    xticklabels({'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'})
    ylabel('Absolute Error')
    
    hold(plot1, 'on')
    hold(plot2, 'on')
    
    % Update delta x
    dx = dx/2;
end

hold(plot1, 'off')
lgd1 = legend(plot1, 'show');
lgd1.Location = 'northwest';

hold(plot2, 'off')
lgd2 = legend(plot2, 'show');
lgd2.Location = 'northwest';
%% Image # 7 for thesis -- Error in Finite Difference Scheme
clear; clc; close all;

a  = 1;     % Amplitude of the sin initial condition
L  = 1;     % Length of the x-interval -- [0, L]
T  = 0.01;  % Length of the t-intercal -- [0, T]
c  = 1/15;  % Diffusion/Viscosity Constant
dx = 1/20;  % This is delta x -- The size of the sub-interval in space
dt = 0.001; % This is delta t -- The size of the sub-interval in time
Col = 4;    % Set the number of columns in the table
ROW = 9;    % This is the number of rows -- this value is treated like a constant * Do not change
numTerms  = 100; % Number of terms in the finite Fourier series
    
% Toggle the chosen initial condition (feel free to define your own!)
u0_1 = @(z) a*sin(pi*z/L);
u0_2 = @(z) 4*z.*(1 - z);

hold on

for j = 1:Col
    
    % Apply the Cole-Hopf Transformation to the initial condition
    x1 = ColeHopfTransformation_Numerical(u0_1, c, L, dx);
    x2 = ColeHopfTransformation_Numerical(u0_2, c, L, dx);

    % Solve the heat equation using finite differences via the Crank-Nicolson scheme
     v1           = HeatEq1D_CrankNicolson(x1, c, L, T, dx, dt);
    [v2, x, t, N] = HeatEq1D_CrankNicolson(x2, c, L, T, dx, dt);

    % Compute the solution to Burgers' equation via the discrete Cole-Hopf transformation
    u1 = BurgersEq1D_Numerical(x1, c, L, T, dx, dt, numTerms);
    u2 = BurgersEq1D_Numerical(x2, c, L, T, dx, dt, numTerms);
    u_approx1 = ColeHopfTransformation_Discrete(v1, c, L, T, dx, dt);
    u_approx2 = ColeHopfTransformation_Discrete(v2, c, L, T, dx, dt);

    % Calculate the approprate abscissae values
    abscissae = N/(ROW + 1) + 1:N/(ROW + 1):N;
    
    % Plot the error
    plot1 = subplot(1, 2, 1);
    txt = ['N = ', num2str(N)];
    plot(abs(u_approx1(abscissae', end) - u1(abscissae', end)), 'DisplayName', txt);
    title(['The Error in Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
    xlabel('x')
    xticks(1:9)
    xticklabels({'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'})
    ylabel('Absolute Error')
    
    plot2 = subplot(1, 2, 2);
    txt = ['N = ', num2str(N)];
    plot(abs(u_approx2(abscissae', end) - u2(abscissae', end)), 'DisplayName', txt);
    title(['The Error in Burgers'' equation for \kappa = ', num2str(c), ' (\itRe = ', num2str(1/c), ')'])
    xlabel('x')
    xticks(1:9)
    xticklabels({'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'})
    ylabel('Absolute Error')
    
    hold(plot1, 'on')
    hold(plot2, 'on')
    
    % Update delta x
    dx = dx/2;
end

hold(plot1, 'off')
lgd1 = legend(plot1, 'show');
lgd1.Location = 'northwest';

hold(plot2, 'off')
lgd2 = legend(plot2, 'show');
lgd2.Location = 'northwest';