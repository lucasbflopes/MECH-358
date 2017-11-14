close all
clc
clear all

% make the mesh
N=50;
xmin=-2;
xmax=2;
ymin=xmin;
ymax=xmax;
x=linspace(xmin,xmax,N);
y=linspace(ymin,ymax,N);
[X, Y]=meshgrid(x,y);

Z=exp(-X.^2-Y.^2);

figure('Position',[30,400,600,400]);
surf(X,Y,Z);
axis tight
set(gca,'FontSize',17);
colormap Hot
% colormap Cool
% colormap Cool works as well
% colormap Gray
title('surf','FontSize',17);
xlabel('x','FontSize',17);
ylabel('y','FontSize',17);

figure('Position',[750,400,600,400]);
contour(X,Y,Z,10,'LineWidth',2);
axis tight
set(gca,'FontSize',17);
% colormap Hot
% colormap Cool
% colormap Cool works as well
colormap jet
title('contour','FontSize',17);
xlabel('x','FontSize',17);
ylabel('y','FontSize',17);
