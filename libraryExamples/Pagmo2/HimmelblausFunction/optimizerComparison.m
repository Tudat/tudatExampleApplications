clc
clear all
close all

saveFolder = '/home/dominic/Software/optimizationBundle/tudatBundle/tudatExampleApplications/libraryExamples/Pagmo2/bin/applications/';


numberOfOptimizers = 16;
numberOfGenerations = 16;

plotInterval = 1;

fitness = cell(numberOfOptimizers,numberOfGenerations);
population = cell(numberOfOptimizers,numberOfGenerations);

for i = 1:numberOfOptimizers
    figure(i)
    counter = 1;
    for j = 1:plotInterval:numberOfGenerations
      fitness{i,j} =  load(strcat(saveFolder,'fitness_',num2str(i-1),'_',num2str(j),'.dat'));
      population{i,j} =  load(strcat(saveFolder,'population_',num2str(i-1),'_',num2str(j),'.dat'));
      subplot(4,4,counter)
      %if( j > 1 )
      scatter(population{i,j}(:,1),population{i,j}(:,2), 10,fitness{i,j},'*')
            %scatter(population{i,j}(:,1)-population{i,j-plotInterval}(:,1),population{i,j}(:,2)-population{i,j-plotInterval}(:,2), 10,fitness{i,j},'*')
      %end
      counter = counter + 1;
      grid on
      %xlim([-5 10])
      %ylim([-5 10])

    end
end