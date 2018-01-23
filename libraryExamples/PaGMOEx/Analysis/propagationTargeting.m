set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

clc
clear all
close all

saveFolder = '/home/dominic/Software/optimizationBundle/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/SimulationOutput/';

load(strcat(saveFolder,'propagationTargetingGridSearch.dat'))
load(strcat(saveFolder,'propagationTargetingGridSearch_x_data.dat'))
load(strcat(saveFolder,'propagationTargetingGridSearch_y_data.dat'))
figure(1)
contour(propagationTargetingGridSearch_x_data,propagationTargetingGridSearch_y_data,propagationTargetingGridSearch')

colorbar

xlabel('Argument of periapsis [deg]')
ylabel('Longitude of asc. node [deg]')

title('Minimum targeting error [m]')
%%
figure(2)

population = cell(8,1);
fitness = cell(8,1);

for k=1:8
    
    %figure(k+1)
    if( k == 1 )
        indexToUse = 1;
    elseif( k == 2 )
        indexToUse = 2;
    elseif( k == 3 )
        indexToUse = 3;
    elseif( k == 4 )
        indexToUse = 4;
    elseif( k == 5 )
        indexToUse = 5;
    elseif( k == 6 )
        indexToUse = 8;
    elseif( k == 7 )
        indexToUse = 12;
    elseif( k == 8 )
        indexToUse = 25;
    end
    
    %figure(k)
    
    j=indexToUse
    population{j} = load(strcat(saveFolder,'population_targetingPropagation_',num2str(j-1),'_',num2str(j-1),'.dat'));
    fitness{j} = load(strcat(saveFolder,'fitness_targetingPropagation_',num2str(j-1),'_',num2str(j-1),'.dat'));
    
    
    subplot(2,4,k)
    contour(propagationTargetingGridSearch_x_data,propagationTargetingGridSearch_y_data,propagationTargetingGridSearch')
    hold on
    scatter(population{indexToUse}(:,1),population{indexToUse}(:,2),25,fitness{indexToUse}(:,1),'*')
    
    xlabel('Argument of periapsis [deg]')
    ylabel('Longitude of asc. node [deg]')
    
    
    min(fitness{indexToUse}(:,1))
    
    title(strcat('Iteration ',{' '},num2str(indexToUse)));
    
end

set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
set(gcf,'PaperPositionMode','auto');

%%

figure(3)

load(strcat(saveFolder,'propagationTargetingGridSearch_pert.dat'))
load(strcat(saveFolder,'propagationTargetingGridSearch_pert_x_data.dat'))
load(strcat(saveFolder,'propagationTargetingGridSearch_pert_y_data.dat'))

contour(propagationTargetingGridSearch_pert_x_data,propagationTargetingGridSearch_pert_y_data,propagationTargetingGridSearch_pert')

colorbar

xlabel('Argument of periapsis [deg]')
ylabel('Longitude of asc. node [deg]')

title('Minimum targeting error [m]')

figure(4)

population_pert = cell(8,1);
fitness_pert = cell(8,1);

for k=1:4
    
    %figure(k+1)
    if( k == 1 )
        indexToUse = 1;
    elseif( k == 2 )
        indexToUse = 2;
    elseif( k == 3 )
        indexToUse = 3;
    elseif( k == 4 )
        indexToUse = 4;
    end
    
    %figure(k)
    
    j=indexToUse
    population_pert{j} = load(strcat(saveFolder,'population_targetingPropagation_pert_',num2str(j-1),'_',num2str(j-1),'.dat'));
    fitness_pert{j} = load(strcat(saveFolder,'fitness_targetingPropagation_pert_',num2str(j-1),'_',num2str(j-1),'.dat'));
    
    min(fitness_pert{j})
    
    subplot(1,4,k)
    contour(propagationTargetingGridSearch_pert_x_data,propagationTargetingGridSearch_pert_y_data,propagationTargetingGridSearch_pert')
    hold on
    scatter(population_pert{indexToUse}(:,1),population_pert{indexToUse}(:,2),25,fitness_pert{indexToUse}(:,1),'*')
    
    xlabel('Argument of periapsis [deg]')
    ylabel('Longitude of asc. node [deg]')
    
    
    min(fitness_pert{indexToUse}(:,1))
    
    title(strcat('Iteration ',{' '},num2str(indexToUse)));
    
end

set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
set(gcf,'PaperPositionMode','auto');
