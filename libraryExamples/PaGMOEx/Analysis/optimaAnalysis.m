set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

clc
close all
clear all

saveFolder = '/home/dominic/Software/optimizationBundle/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/SimulationOutput/';

optima_2_dim = cell(3,1);
plotOptions = cell(3,1);
plotOptions{1} = 'b';
plotOptions{2} = 'r';
plotOptions{3} = 'k';

for k=1:3
    optima_2_dim{k} = load(strcat(saveFolder,'cec2013Optima_0_',num2str(k-1),'.dat'));
end

figure(1)
for k=1:3
    [Y,I]=min(optima_2_dim{k}');
    %     for i=1:11
    %         numberOfOptima(i) = sum(I==i);
    for j=1:28
        subplot(4,7,j)
        if( j == 1)
            plotHandle{k}=plot(optima_2_dim{k}(j,:),strcat(plotOptions{k},'-*'));
        else
            plot(optima_2_dim{k}(j,:),strcat(plotOptions{k},'-*'))
        end
        hold on
        grid on
        scatter(I(j),Y(j),50*k,strcat(plotOptions{k},'o'))
        xlim([0 12])
        title(strcat('Problem ',{' '},num2str(j)))
        if(rem(j,7)==1)
            ylabel('Minimum [-]')
        end
        if(j>21)
            xlabel('Optimizer [-]')
        end
        
    end
end

subplot(4,7,1)
legend([plotHandle{1} plotHandle{2} plotHandle{3}],'Pop: 16, Gen: 64', 'Pop: 32, Gen: 32', 'Pop: 64, Gen: 16'  )

set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 60 40]);
set(gcf,'PaperPositionMode','auto');

suptitle('Optima, CEC2013 2-D test problems')

%%
for k=1:3
    optima_5_dim{k} = load(strcat(saveFolder,'cec2013Optima_1_',num2str(k-1),'.dat'));
end

figure(2)
for k=1:3
    [Y,I]=min(optima_5_dim{k}');
    for j=1:28
        subplot(4,7,j)
        if( j == 1)
            plotHandle{k}=plot(optima_5_dim{k}(j,:),strcat(plotOptions{k},'-*'));
        else
            plot(optima_5_dim{k}(j,:),strcat(plotOptions{k},'-*'))
        end
        
        hold on
        grid on
        scatter(I(j),Y(j),50*k,strcat(plotOptions{k},'o'))
        xlim([0 12])
        title(strcat('Problem ',{' '},num2str(j)))
        if(rem(j,7)==1)
            ylabel('Minimum [-]')
        end
        if(j>21)
            xlabel('Optimizer [-]')
        end
        
    end
end


subplot(4,7,1)
legend([plotHandle{1} plotHandle{2} plotHandle{3}],'Pop: 16, Gen: 64', 'Pop: 32, Gen: 32', 'Pop: 64, Gen: 16' )

suptitle('Optima, CEC2013 5-D test problems')

set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 60 40]);
set(gcf,'PaperPositionMode','auto');


%%
for k=1:3
    optima_10_dim{k} = load(strcat(saveFolder,'cec2013Optima_2_',num2str(k-1),'.dat'));
end

figure(3)
for k=1:3
    [Y,I]=min(optima_10_dim{k}');
    for j=1:28
        subplot(4,7,j)
        if( j == 3)
            plotHandle{k}=plot(optima_10_dim{k}(j,:),strcat(plotOptions{k},'-*'));
        else
            plot(optima_10_dim{k}(j,:),strcat(plotOptions{k},'-*'))
        end
        hold on
        grid on
        scatter(I(j),Y(j),50*k,strcat(plotOptions{k},'o'))
        xlim([0 12])
        title(strcat('Problem ',{' '},num2str(j)))
        if(rem(j,7)==1)
            ylabel('Minimum [-]')
        end
        if(j>21)
            xlabel('Optimizer [-]')
        end
        
    end
end

subplot(4,7,1)
legend([plotHandle{1} plotHandle{2} plotHandle{3}],'Pop: 16, Gen: 64', 'Pop: 32, Gen: 32', 'Pop: 64, Gen: 16'  )

suptitle('Optima, CEC2013 10-D test problems')

set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 60 40]);
set(gcf,'PaperPositionMode','auto');

