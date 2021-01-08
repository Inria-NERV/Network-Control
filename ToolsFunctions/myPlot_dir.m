function h = myPlot_dir (vX, vY , mySymbol)
myColor = 'r';

%function [h1, h2] = plot_dir (vX, vY)
%Plotting x-y variables with direction indicating vector to the next element.
%Example
%   vX = linspace(0,2*pi, 10)';
%   vY = sin (vX);
%   plot_dir(vX, vY);

rMag = 0.5;

% Length of vector
lenTime = length(vX);

% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;

% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);

% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);

% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;

% make plot 
%h1 = plot (vX, vY, 'Color' , myColor , 'LineStyle' , '-'); hold on;
if strcmp(myColor , 'r')
    colormap jet
elseif strcmp(myColor , 'b')
    colormap jet
end
z = 1:lenTime;
h = scatter (vX, vY, 17 , z, 'filled' , mySymbol); hold on;
% add arrows 
% h2 = quiver (vXQ0,vYQ0, vPx, vPy, 0, myColor ,'LineWidth' , 2 );%'AutoScale' , 'off' , 'AutoScaleFactor' , 5);
grid on; hold off
axis equal
