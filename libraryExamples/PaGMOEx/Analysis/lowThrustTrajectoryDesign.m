%
% This script processes the results of the Earth-Mars Lambert targeting
% multi-objective  optimization, run by the
% multiObjectiveEarthMarsTransferExample.cpp Tudat/Pagmo2 example.
%

set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

clc
clear all
close all

% Load grid search data
% saveFolder = '../SimulationOutput/';
hodographicShaping = importdata('../SimulationOutput/hodographicShapingGridSearch.dat');
hodographicShapingHigherOrder = importdata('../SimulationOutput/hodographicShapingResultsHigherOrder.dat');
hodographicShapingOneRevolution = importdata('../SimulationOutput/hodographicShapingOneRevolution.dat');

for i = 1:401 %size(hodographicShaping,1)
    departureDate(i,1) = hodographicShaping(i,3);
end
for i = 1:301 %size(hodographicShaping,1)
    TOF(i,1) = hodographicShaping( ( i - 1 ) * 401 + 1,2);
end
for i = 1:301
    for j = 1:401
        deltaV(i,j) = hodographicShaping( ( i - 1 ) * 401 + j, 4 );
    end
end

for i = 1:6 %size(hodographicShaping,1)
    departureDateOneRevolution(i,1) = hodographicShapingOneRevolution(i,3);
end
for i = 1:21 %size(hodographicShaping,1)
    TOFoneRevolution(i,1) = hodographicShapingOneRevolution( ( i - 1 ) * 6 + 1,2);
end
for i = 1:21
    for j = 1:6
        deltaVoneRevolution(i,j) = hodographicShapingOneRevolution( ( i - 1 ) * 6 + j, 4 );
    end
end


figure(1);
hold on;
contourf( departureDate', TOF', deltaV / 1e3, 3000, 'LineColor', 'None' );
% plot( hodographicShapingOneRevolution(:,3), hodographicShapingOneRevolution(:,2), '.k', 'markerSize', 10);
grid;
hold off;
ylabel('TOF [days]', 'fontSize', 11, 'Interpreter', 'Latex');
xlabel('Departure date [MJD2000]', 'fontSize', 11, 'Interpreter', 'Latex');

c1=zeros(45,3); 
c1(1:7,1)= 2/256; 
c1(1:7,2)= 29/256;
c1(1:7,3)= 242/256;
c1(8,1)= 0/256;
c1(8,2)= 136/256;
c1(8,3)= 255/256; 
c1(9:10,1)= 1/256;
c1(9:10,2)= 171/256;
c1(9:10,3)= 255/256; 
c1(11:15,1)=1/256;  
c1(11:15,2)=213/256;
c1(11:15,3)=253/256;
c1(16:20,1)=55/256; 
c1(16:20,2)=255/256;
c1(16:20,3)=200/256;
c1(21:40,1)=158/256;
c1(21:40,2)=255/256;
c1(21:40,3)=98/256;
c1(41:45,1)= 203/256;
c1(41:45,2)= 2/256;
c1(41:45,3)= 0/256;
colormap(c1);
c = colorbar;
caxis([0 45]);
set(c,'YTick',[7 8 10 15 20 40 45])
title(c,'$\Delta V$ [km/s]', 'Interpreter', 'Latex', 'fontSize', 10);



%% High-order solution

for i = 1:6 %size(hodographicShaping,1)
    departureDateHighOrder(i,1) = hodographicShapingHigherOrder(i,3);
end
for i = 1:21 %size(hodographicShaping,1)
    TOFhighOrder(i,1) = hodographicShapingHigherOrder( ( i - 1 ) * 6 + 1,2);
end
for i = 1:21
    for j = 1:6
        deltaVhighOrder(i,j) = hodographicShapingHigherOrder( ( i - 1 ) * 6 + j, 4 );
    end
end


figure(2);
subplot(1,2,1);
hold on;
contourf( departureDateHighOrder', TOFhighOrder', deltaVhighOrder / 1e3, 3000, 'LineColor', 'None' );
plot( hodographicShapingHigherOrder(:,3), hodographicShapingHigherOrder(:,2), '.k', 'markerSize', 10);
grid;
hold off;
ylabel('TOF [days]', 'fontSize', 11, 'Interpreter', 'Latex');
xlabel('Departure date [MJD2000]', 'fontSize', 11, 'Interpreter', 'Latex');

colormap(c1);
c = colorbar;
caxis([0 45]);
% xlim([departureDateHighOrder(1,1) 7350]);
set(c,'YTick',[7 8 10 15 20 40 45])
title(c,'$\Delta V$ [km/s]', 'Interpreter', 'Latex', 'fontSize', 10);
title('High-order solution');

figure(2);
subplot(1,2,2);
hold on;
contourf( departureDateOneRevolution', TOFoneRevolution', deltaVoneRevolution / 1e3, 3000, 'LineColor', 'None' );
plot( hodographicShapingOneRevolution(:,3), hodographicShapingOneRevolution(:,2), '.k', 'markerSize', 10);
grid;
hold off;
ylabel('TOF [days]', 'fontSize', 11, 'Interpreter', 'Latex');
xlabel('Departure date [MJD2000]', 'fontSize', 11, 'Interpreter', 'Latex');
title('Low-order solution');

colormap(c1);
c = colorbar;
caxis([0 45]);
set(c,'YTick',[7 8 10 15 20 40 45])
xlim([departureDateHighOrder(1) departureDateHighOrder(end)]);
ylim([TOFhighOrder(1) TOFhighOrder(end)])
title(c,'$\Delta V$ [km/s]', 'Interpreter', 'Latex', 'fontSize', 10);





% c1(51:99,1)=1; %yellow
% c1(51:99,2)=1; %yellow
% colormap(c)
% colorbar

% 
% % Create grid search plot
% contour(himmelBlauGridSearch_x_data,himmelBlauGridSearch_y_data,himmelBlauGridSearch',50)
% xlabel('x [-]')
% ylabel('y [-]')
% title('Himmelblau function value [-]')
% colorbar
% 
% set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
% set(gcf,'PaperPositionMode','auto');
% 
% pause(0.1)
% saveas(figure(1),'himmelblauGridSearch','png');
% %%
% % Plot data for 6 generations
% figure(2)
% population = cell(1,1);
% fitness = cell(1,1);
% for k=1:6
%     
%     % Specify current generation
%     if( k == 1 )
%         indexToUse = 1;
%     elseif( k == 2 )
%         indexToUse = 5;
%     elseif( k == 3 )
%         indexToUse = 10;
%     elseif( k == 4 )
%         indexToUse = 20;
%     elseif( k == 5 )
%         indexToUse = 50;
%     elseif( k == 6 )
%         indexToUse = 100;
%     end
%     
%     % Load data for current generation, for each of the three optimizers
%     for j=indexToUse
%         population{j} = load(strcat(saveFolder,'population_himmelblau_',num2str(j),'.dat'));
%         fitness{j} = load(strcat(saveFolder,'fitness_himmelblau_',num2str(j),'.dat'));
%     end
%     
%     % Plot population for current generation/optimizer on top of porkchop
%     subplot(2,3,k)
%     contour(himmelBlauGridSearch_x_data,himmelBlauGridSearch_y_data,himmelBlauGridSearch',50)
%     xlabel('x [-]')
%     ylabel('y [-]')
%     hold on
%     scatter(population{indexToUse}(:,1),population{indexToUse}(:,2),25,fitness{indexToUse}(:,1),'*')
%     title(strcat('Iteration ',{' '},num2str(indexToUse),{' '},', Minimum=',num2str(min(fitness{j}(:,1)),4)));
%     
% end
% %%
% set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
% set(gcf,'PaperPositionMode','auto');
% 
% pause(0.1)
% saveas(gcf,strcat('himmelblauOptimization'),'png');
