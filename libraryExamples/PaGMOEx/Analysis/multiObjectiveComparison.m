set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

clc
clear all
close all

saveFolder = '/home/dominic/Software/optimizationBundle/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/SimulationOutput/';


numberOfOptimizers = 1;
numberOfGenerations = 32;

plotInterval = 2;

fitness = cell(numberOfOptimizers,numberOfGenerations);
population = cell(numberOfOptimizers,numberOfGenerations);

counter = 1;
i = 1;
for j = 1:plotInterval:numberOfGenerations
    fitness{i,j} =  load(strcat(saveFolder,'fitness_mo_mga_EVEEJ_',num2str(j),'.dat'));
    population{i,j} =  load(strcat(saveFolder,'population_mo_mga_EVEEJ_',num2str(j),'.dat'));
    
    figure(1)
    subplot(4,4,counter)
    scatter(fitness{i,j}(:,1),fitness{i,j}(:,2), 10,'*')
    
    if(rem(counter,4)==1)
        ylabel('Travel time [days]')
    end
    if(counter > 12 )
        xlabel('\Delta V [km/s]')
    end
    title(strcat('Iteration ',{' '},num2str(j)));
    
    grid on
    for k=1:size(population{i,j},2)
        figure(k+1)
        subplot(4,4,counter)
        if( k == 1 )
            scatter(fitness{i,j}(:,1),fitness{i,j}(:,2), 10,(population{i,j}(:,k)-2451545)/365,'*')
        else
            scatter(fitness{i,j}(:,1),fitness{i,j}(:,2), 10,population{i,j}(:,k),'*')
        end
        
        if(rem(counter,4)==1)
            ylabel('Travel time [days]')
        end
        if(counter > 12 )
            xlabel('Delta V [km/s]')
        end
        
        title(strcat('Iteration ',num2str(j)));
        
        colorbar
        grid on
    end
    counter = counter + 1;
    grid on    
end

for k=1:size(population{1,1},2)+1
    
    set(figure(k), 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
    set(figure(k),'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
    set(figure(k),'PaperPositionMode','auto');
end

for k=1:size(population{i,j},2)
    figure(k+1)
    if( k== 1)
        suptitle('Color scale: departure date [years since J2000]')
    elseif( k== 2)
        suptitle('Color scale: leg 1 (E-V) travel time [days]')
    elseif( k== 3)
        suptitle('Color scale: leg 2 (V-E) travel time [days]')
    elseif( k== 4)
        suptitle('Color scale: leg 3 (E-E) travel time [days]')
    elseif( k== 5)
        suptitle('Color scale: leg 4 (E-J) travel time [days]')
    end
end