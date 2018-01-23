set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

clc
clear all
close all

saveFolder = '/home/dominic/Software/optimizationBundle/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/SimulationOutput/';

load(strcat(saveFolder,'porkchopEarthMars.dat'))
load(strcat(saveFolder,'porkchopEarthMars_x_data.dat'))
load(strcat(saveFolder,'porkchopEarthMars_y_data.dat'))

contour((porkchopEarthMars_x_data-2451545)/365,porkchopEarthMars_y_data,porkchopEarthMars')
xlabel('Departure date [years since J2000]')
ylabel('Travel time [days]')
title('Porkchop plot, two-shot impulsive Earth-Mars transfer [m/s]')

set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
set(gcf,'PaperPositionMode','auto');
%%
optimizers = cell(3,1);
optimizers{1} = 'NSGA2';
optimizers{2} = 'MOEAD';
optimizers{3} = 'IHS';

for k=1:6
    population = cell(11,16);
    fitness = cell(11,16);
    figure(k+1)
    if( k == 1 )
        indexToUse = 1;
    elseif( k == 2 )
        indexToUse = 5;
    elseif( k == 3 )
        indexToUse = 10;
    elseif( k == 4 )
        indexToUse = 20;
    elseif( k == 5 )
        indexToUse = 50;
    elseif( k == 6 )
        indexToUse = 100;
    end
    
    for i=1:3
        for j=indexToUse
            population{i,j} = load(strcat(saveFolder,'population_mo_EarthMars_',num2str(j-1),'_',num2str(i-1),'.dat'));
            fitness{i,j} = load(strcat(saveFolder,'fitness_mo_EarthMars_',num2str(j-1),'_',num2str(i-1),'.dat'));
        end
        subplot(2,3,i)
        contour((porkchopEarthMars_x_data-2451545)/365,porkchopEarthMars_y_data,porkchopEarthMars')
        hold on
        scatter((population{i,indexToUse}(:,1)-2451545)/365,population{i,indexToUse}(:,2),25,fitness{i,indexToUse}(:,1),'*')
        title(optimizers{i})
        xlabel('Departure date [years since J2000]')
        ylabel('Travel time [days]')
        subplot(2,3,i+3)
        scatter(fitness{i,indexToUse}(:,1),fitness{i,indexToUse}(:,2),'*')
        grid on
        xlabel('Delta V [km/s]','interpreter','latex')
        ylabel('Travel time [days]')
        title(strcat('Minimum Delta V:',{' '},num2str(min(fitness{i,indexToUse}(:,1)))),'interpreter','latex')
        
    end
    suptitle(strcat('Iteration ',{' '},num2str(indexToUse)));
    
    set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
    set(gcf,'PaperPositionMode','auto');
end
