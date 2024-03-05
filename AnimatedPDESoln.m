function AnimatedPDESoln(u, u0, c, x, t, L, a, filename)
% This function animatedPDESoln, animates the solution to the given PDE. No objects get passed 
% back to the user by this method, rather it simply creates the animation, closes it when it is finished 
% running (i.e. max time is reached -- all of the vector t has been traversed), and finally saves a gif of
% the solution. This function needs to be tailored to the chosen PDE.
%
%                   u = The solution to the given PDE          -- matrix of real numbers
%                  x0 = The initial condition passed as vector -- positive real numbers
%                   c = The diffusion/viscosity constant       -- positive real number
%                   x = The x-values that the PDE is solved at -- vector of real numbers
%                   t = The t-values that the PDE is solved at -- vector of real numbers
%                   L = The end point of the interval in space -- [0, L] (positive real number)
%                   A = The amplitude of the sin initial condition (if chosen) -- real number
%            filename = The chosen filename passed to the function as a string
%
% @author Jarren Ralf

% Set the initial condition, instantiate a plot object, and then set the graph axis
  y = u0(x);
fig = figure;
  h = plot(x, y);
axis([0 L 0 a]);

% Set the X-Data and Y-Data source for the plot object to be x and y respectively
h.XDataSource = 'x';
h.YDataSource = 'y';

% Retreive the number of columns of the solution matrix
numColumns = size(u, 2); 

% Initialize the cell array that will store frame structures corresponding to the figure
images = cell(1, numColumns);

% Animate and create a gif!
for i = 1:numColumns
    
    % Pass the updated function values to y, i.e. the solution changing in time
    y = u(:, i); 
    
    % Updates the data source, i.e. the y-values from above and label the axis
    plot(x, y);
    axis([0 L 0 a]);
    refreshdata
    title({['The Solution of Burgers'' equation for \nu = ', num2str(c), ' (\itRe = ', num2str(1/c), ...
            ')']; ['t = ', num2str(t(i), '%5.3f')]})
    xlabel('x')
    ylabel('u(x, t)')
    
    % Draws the updated graphics on the figure immediately, then store the frames in the cell array
    drawnow
    frame = getframe(fig);
    images{i} = frame2im(frame);
end

% Close the figure 
close;

% Create the gif and save as the chosen filename in the current directory
for i = 1:numColumns
    [A, map] = rgb2ind(images{i}, 256);
    if i == 1
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end