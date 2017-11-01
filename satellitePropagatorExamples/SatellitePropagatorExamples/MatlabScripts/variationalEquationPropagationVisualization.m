clc
close all
clear all

%%
%%%%% DEFINE DATA DIRECTORY AND LOAD FILES
dataDirectory = '../SimulationOutput/PerturbedSatelliteVariationalExample/';

sensitivityMatrix = load(strcat(dataDirectory,'singlePerturbedSatelliteSensitivityHistory.dat'));
stateTransitionMatrix = load(strcat(dataDirectory,'singlePerturbedSatelliteStateTransitionHistory.dat'));

%%

time = (stateTransitionMatrix(:,1)-stateTransitionMatrix(1,1))/3600;

colors = cell(3,1);
colors{1} = 'b';
colors{2} = 'r';
colors{3} = 'k';
colors{4} = 'g';
colors{5} = 'c';

styles = cell(3,1);
styles{1} = '-';
styles{2} = '--';
styles{3} = ':';

lineWidths = cell(3,1);
lineWidths{1} = 1;
lineWidths{2} = 1;
lineWidths{3} = 2;

figure(1)
for i=1:2
   subplot(2,2,i)
   for j=1:3
       for k=1:3
           h(j,k)=plot(time,stateTransitionMatrix(:,1+3*(i-1)+(k-1)*6+j),strcat(colors{j},styles{k}),'LineWidth',lineWidths{k})
           hold on
       end
   end
   
   if( i == 1 )
      legend([h(1,1), h(2,1),h(3,1)],{'$x_{0}$','$y_{0}$','$z_{0}$'},'interpreter','latex','Location','NorthWest')
      xlabel('t [h]','interpreter','latex')
      ylabel('$\partial \mathbf{r}/\partial\mathbf{r}_{0}$ [m/m]','interpreter','latex')
   elseif( i == 2 )
      legend([h(1,1), h(1,2),h(1,3)],{'$x$','$y$','$z$'},'interpreter','latex','Location','NorthWest')
         xlabel('t [h]','interpreter','latex')
      ylabel('$\partial \mathbf{r}/\partial\mathbf{v}_{0}$ [m/(m/s)]','interpreter','latex')
   end
   grid on
   
   subplot(2,2,i + 2)
   for j=1:3
       for k=1:3
           plot(time,stateTransitionMatrix(:,19+3*(i-1)+(k-1)*6+j),strcat(colors{j},styles{k}),'LineWidth',lineWidths{k})
           hold on
       end
   end
   
   if( i == 1 )
       legend([h(1,1), h(2,1),h(3,1)],{'$x_{0}$','$y_{0}$','$z_{0}$'},'interpreter','latex','Location','NorthEast')
       xlabel('t [h]','interpreter','latex')
       ylabel('$\partial \mathbf{v}/\partial\mathbf{r}_{0}$ [(m/s)/m]','interpreter','latex')
   elseif( i == 2 )
       legend([h(1,1), h(1,2),h(1,3)],{'$x$','$y$','$z$'},'interpreter','latex','Location','NorthEast')
       xlabel('t [h]','interpreter','latex')
       ylabel('$\partial \mathbf{v}/\partial\mathbf{v}_{0}$ [(m/s)/(m/s)]','interpreter','latex')
   end
   grid on
end
suptitle('State transition matrix, plotted per block')
%%
close all
figure(2)
   subplot(2,2,1)
   for j=1:2
       for k=1:3
           h(j,k)=plot(time,sensitivityMatrix(:,1+(k-1)*7+j),strcat(colors{j},styles{k}),'LineWidth',lineWidths{k});
           hold on
       end
   end
   grid on
   
   subplot(2,2,3)
   for j=1:2
       for k=1:3
           h(j,k)=plot(time,sensitivityMatrix(:,22+(k-1)*7+j),strcat(colors{j},styles{k}),'LineWidth',lineWidths{k});
           hold on
       end
   end
   grid on
   
   
   subplot(2,2,2)
   for j=1:5
       for k=1:3
           h(j,k)=plot(time,sensitivityMatrix(:,3+(k-1)*7+j),strcat(colors{j},styles{k}),'LineWidth',lineWidths{k});
           hold on
       end
   end
   grid on
   
   subplot(2,2,4)
   for j=1:5
       for k=1:3
           h(j,k)=plot(time,sensitivityMatrix(:,24+(k-1)*7+j),strcat(colors{j},styles{k}),'LineWidth',lineWidths{k});
           hold on
       end
   end
   grid on


%%
for j=1:2
    set( figure(j), 'Units', 'normalized', 'Position', [0,0,1,1]);
    set( figure(j),'PaperUnits','centimeters','PaperPosition',[0 0 60 40]);
end