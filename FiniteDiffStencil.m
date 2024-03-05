%% 2-Dimensional Stencil (PDE)
clear; clc; close all;

% horizontal line
x1 = linspace(-1, 1);
y1 = 0*x1;

% vertical line
y2 = linspace(0, .5);
x2 = 0*y2;

subplot(1, 2, 1)
hold on
    axis([-1.1 1.1 -0.4 .9])
    set(gca,'visible','off')
    plot(x1, y1,  'k');
    plot(x2, y2,  'k');
    plot(-1,  0, 'ok', 'MarkerFaceColor', 'r');
    text(-1,  -0.01, 'i-1,j', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot( 0,  0, 'ok', 'MarkerFaceColor', 'r');
    text( 0,  -0.01, 'i,j', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot( 1,  0, 'ok', 'MarkerFaceColor', 'r');
    text( 1,  -0.01, 'i+1,j', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot( 0,  .5, 'ok', 'MarkerFaceColor', 'r');
    text( 0,  .52, 'i,j+1', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
hold off

subplot(1, 2, 2)
hold on
    axis([-1.1 1.1 -0.4 .9])
    set(gca,'visible','off')
    plot(x1, y1,  'k');
    plot(x2, y2,  'k');
    plot(x1, y1 + .5,  'k');
    plot(-1,  0, 'ok', 'MarkerFaceColor', 'r');
    text(-1,  -0.01, 'i-1,j', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot( 0,  0, 'ok', 'MarkerFaceColor', 'r');
    text( 0,  -0.01, 'i,j', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot( 1,  0, 'ok', 'MarkerFaceColor', 'r');
    text( 1,  -0.01, 'i+1,j', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot( 0,  .5, 'ok', 'MarkerFaceColor', 'r');
    text( 0,  .52, 'i,j+1', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    plot(-1,  .5, 'ok', 'MarkerFaceColor', 'r');
    text(-1,  .52, 'i-1,j+1', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    plot( 1,  .5, 'ok', 'MarkerFaceColor', 'r');
    text( 1,  .52, 'i+1,j+1', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
hold off
%% 5-Point Stencil
clear; clc; close all;

% horizontal line
x1 = linspace(-1, 1);
y1 = 0*x1;

hold on
    axis([-1.1 1.1 -0.4 .9])
    set(gca,'visible','off')
    plot(x1, y1,  'k');
    plot(-1  ,  0   , 'ok', 'MarkerFaceColor', 'r');
    text(-1  , -0.01, 'x_0 - 2h', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot(-0.5,  0   , 'ok', 'MarkerFaceColor', 'r');
    text(-0.5, -0.01, 'x_0 - h', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot( 0  ,  0   , 'ok', 'MarkerFaceColor', 'r');
    text( 0  , -0.01, 'x_0', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot( 0.5,  0   , 'ok', 'MarkerFaceColor', 'r');
    text( 0.5, -0.01, 'x_0 + h', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot( 1  ,  0   , 'ok', 'MarkerFaceColor', 'r');
    text( 1  , -0.01, 'x_0 + 2h', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off