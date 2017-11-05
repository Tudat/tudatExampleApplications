clc
close all
clear all

% Define folder in which the data files can be found
dataDirectory = '../SimulationOutput/';

%%
%%%%% PLOT TWO-SATELLITE PROPAGATION 

figure(1)

load(strcat(dataDirectory,'TwoSatelliteExample/asterixPropagationHistory.dat'));
load(strcat(dataDirectory,'TwoSatelliteExample/obelixPropagationHistory.dat'));
scatter3(asterixPropagationHistory(:,2),asterixPropagationHistory(:,3),asterixPropagationHistory(:,4))
hold on
scatter3(obelixPropagationHistory(:,2),obelixPropagationHistory(:,3),obelixPropagationHistory(:,4))
axis equal

legend('Asterix orbit','Obelix orbit')

xlabel('x_{J2000} [m]')
ylabel('y_{J2000} [m]')
zlabel('z_{J2000} [m]')
%%
%%%%% PLOT GALILEO CONSTELLATION PROPAGATION

figure(2)

for i=1:30
    currentGalileo = load(strcat(dataDirectory,'GalileoConstellationExample/galileoSatellite',num2str(i),'.dat'));
    scatter3(currentGalileo(:,2),currentGalileo(:,3),currentGalileo(:,4))
    hold on
    axis equal
end
xlabel('x_{J2000} [m]')
ylabel('y_{J2000} [m]')
zlabel('z_{J2000} [m]')

%%
%%%%% PLOT APOLLO CAPSULE ENTRY BEHAVIOUR
figure(3)

load(strcat(dataDirectory,'ApolloCapsuleExample/apolloPropagationHistory.dat'));
load(strcat(dataDirectory,'ApolloCapsuleExample/apolloDependentVariableHistory.dat'));


subplot(2,3,1)
plot(apolloDependentVariableHistory(:,1),apolloDependentVariableHistory(:,3));
grid on
xlabel('Time since propagation start [s]')
ylabel('Earth altitude [km]')

subplot(2,3,2)
speed = sqrt(apolloPropagationHistory(:,5).^2+apolloPropagationHistory(:,6).^2+apolloPropagationHistory(:,7).^2);
plot(apolloPropagationHistory(:,1),speed);
grid on
xlabel('Time since propagation start [s]')
ylabel('Earth-centered inertial speed [m/s]')

subplot(2,3,3)
plot(apolloDependentVariableHistory(:,1),apolloDependentVariableHistory(:,2));
grid on
xlabel('Time since propagation start [s]')
ylabel('Mach number [-]')


subplot(2,3,4)
semilogy(apolloDependentVariableHistory(:,1),apolloDependentVariableHistory(:,4));
grid on
xlabel('Time since propagation start [s]')
ylabel('Aerodynamic acceleration [m/s^2]')

subplot(2,3,5)
plot(apolloDependentVariableHistory(:,1),-apolloDependentVariableHistory(:,5));
grid on
xlabel('Time since propagation start [s]')
ylabel('Drag coefficient [-]')

subplot(2,3,6)
plot(apolloDependentVariableHistory(:,1),-apolloDependentVariableHistory(:,7));
grid on
xlabel('Time since propagation start [s]')
ylabel('Lift coefficient [-]')


%%
%%%% PLOT SOLAR SYSTEM DYNAMICS SIMULATION

propagatedBodies = cell(6,1);

propagatedBodies{1} = 'Moon';
propagatedBodies{2} = 'Earth';
propagatedBodies{3} = 'Mars';
propagatedBodies{4} = 'Venus';
propagatedBodies{5} = 'Mercury';
propagatedBodies{6} = 'Sun';

centralBodies{1} = 'Earth';
centralBodies{2} = 'Sun';
centralBodies{3} = 'Sun';
centralBodies{4} = 'Sun';
centralBodies{5} = 'Sun';
centralBodies{6} = 'Barycenter';

for i=1:6
    barycentricSolarSystemOrbits{i} = load(strcat(dataDirectory,'InnerSolarSystemPropagationExample/innerSolarSystemPropagationHistory',propagatedBodies{i},'0.dat'));
    figure(4)
    scatter3(barycentricSolarSystemOrbits{i}(:,2),barycentricSolarSystemOrbits{i}(:,3),barycentricSolarSystemOrbits{i}(:,4))
    hold on
    
    figure(5)
    subplot(3,2,i)
    scatter3(barycentricSolarSystemOrbits{i}(:,2),barycentricSolarSystemOrbits{i}(:,3),barycentricSolarSystemOrbits{i}(:,4))
    title(strcat('Barycentric orbit of ',propagatedBodies{i}));
    
    localBarycentricSolarSystemOrbits{i} = load(strcat(dataDirectory,'InnerSolarSystemPropagationExample/innerSolarSystemPropagationHistory',propagatedBodies{i},'1.dat'));
    figure(6)
    subplot(3,2,i)
    scatter3( localBarycentricSolarSystemOrbits{i}(:,2), localBarycentricSolarSystemOrbits{i}(:,3), localBarycentricSolarSystemOrbits{i}(:,4))
    title(strcat('Orbit of ',propagatedBodies{i},' w.r.t.',centralBodies{i}));
    
    
end

figure(4)
axis equal
legend('Moon','Earth','Mars','Venus','Mercury','Sun')
title('Inner solar system')
xlabel('x_{ECLIPJ2000} [m]')
ylabel('y_{ECLIPJ2000} [m]')
zlabel('z_{ECLIPJ2000} [m]')

figure(5)
for i=1:6
    subplot(3,2,i)
    
    xlabel('x_{ECLIPJ2000} [m]')
    ylabel('y_{ECLIPJ2000} [m]')
    zlabel('z_{ECLIPJ2000} [m]')
end

figure(6)
for i=1:6
    subplot(3,2,i)
    
    xlabel('x_{ECLIPJ2000} [m]')
    ylabel('y_{ECLIPJ2000} [m]')
    zlabel('z_{ECLIPJ2000} [m]')
end

%%
%%%%% PLOT ORBIT OF (UN)PERTURBED EARTH ORBITER, AND DIFFERENCE BETWEEN THE TWO

figure(7)
perturbedSatellite = load(strcat(dataDirectory,'PerturbedSatelliteExample/singlePerturbedSatellitePropagationHistory.dat'));
unperturbedSatellite = load(strcat(dataDirectory,'UnperturbedSatelliteExample/singleSatellitePropagationHistory.dat'));

for i=1:3
    
    subplot(3,2,2*(i-1)+1)
    plot(perturbedSatellite(:,1),perturbedSatellite(:,i+1))
    hold on
    grid on
    
    plot(unperturbedSatellite(:,1),unperturbedSatellite(:,i+1),'r--')
    subplot(3,2,2*i)
    plot(perturbedSatellite(:,1),perturbedSatellite(:,i+1)-unperturbedSatellite(:,i+1))
    grid on

    if( i == 3 )
        subplot(3,2,2*(i-1)+1)
        xlabel('Time since epoch [s]')
        
        subplot(3,2,2*i)
        xlabel('Time since epoch [s]')
        
    end
    
    if( i == 1 )
        subplot(3,2,2*(i-1)+1)
        title('Absolute orbit')
        legend('Perturbed','Unperturbed')
        
        subplot(3,2,2*i)
        title('Orbit difference')
    end
    subplot(3,2,2*(i-1)+1)
    if( i == 1 )
        ylabel('x-component [m]')
    elseif( i == 2 )
        ylabel('y-component [m]')
        
    else
        ylabel('z-component [m]')
        
    end
    
end
%%
%%%%% PLOT RESULTS OF THRUST EXAMPLE (THRUST ALONG VELOCITY VECTOR)

figure(8)
load(strcat(dataDirectory,'ThrustAlongVelocityExample/velocityVectorThrustExample.dat'));

subplot(1,2,1)
plot(velocityVectorThrustExample(:,2),velocityVectorThrustExample(:,3))
axis equal
grid on
title('Cartesian position w.r.t. Earth')
xlabel('x-component [m]')
ylabel('y-component [m]')

subplot(1,2,2)
plot(velocityVectorThrustExample(:,1)-velocityVectorThrustExample(1,1),velocityVectorThrustExample(:,8))
grid on
xlabel('Time since propagation start [s]')
ylabel('Vehicle mass [kg]')
xlim([0 14*86400])

%%
%%%%% PLOT RESULTS OF THRUST EXAMPLE (THRUST FROM FILE)

figure(9)
load(strcat(dataDirectory,'ThrustFromFileExample/thrustExampleFromFilePropagationHistory.dat'));

plot3(thrustExampleFromFilePropagationHistory(:,2),thrustExampleFromFilePropagationHistory(:,3),thrustExampleFromFilePropagationHistory(:,4))
axis equal
grid on
title('Cartesian position w.r.t. Earth')
xlabel('x-component [m]')
ylabel('y-component [m]')
zlabel('z-component [m]')

load('thrustExampleFromFileDependentVariableHistory.dat');

colors=cell(3,1);
colors{1}='r';
colors{2}='b';
colors{3}='k';

figure(10)
for i=1:3
   subplot(1,2,1)
   plot(thrustExampleFromFileDependentVariableHistory(:,1),thrustExampleFromFileDependentVariableHistory(:,i+1),colors{i},'LineWidth',2)
   hold on
   grid on
   
   subplot(1,2,2)
   plot(thrustExampleFromFileDependentVariableHistory(:,1),thrustExampleFromFileDependentVariableHistory(:,i+13),colors{i},'LineWidth',2)
   hold on
   grid on
end
subplot(1,2,1)
xlabel('Time [s]')
ylabel('Thrust acceleration (inertial) [m/s^2]')

subplot(1,2,2)
xlabel('Time [s]')
xlabel('Time [s]')
ylabel('Thrust force (LVLH) [N]')
legend('x-component','y-component','z-component','Location','NorthWest')



%%
%%%% SAVE ALL FIGURES TO FILES

for j=1:10
    set( figure(j), 'Units', 'normalized', 'Position', [0,0,1,1]);
    set( figure(j),'PaperUnits','centimeters','PaperPosition',[0 0 60 40]);
end

pause(2.0)

for j=1:10
    figure(j)    
    if(j==1)
        saveas(gcf,strcat('asterixObelixOrbits'),'png');
    elseif(j==2)
        saveas(gcf,strcat('galileoConstellationOrbits'),'png');
    elseif(j==3)
        suptitle('Capsule entry properties')
        saveas(gcf,strcat('capsuleEntryProperties'),'png');
    elseif(j==4)
        saveas(gcf,strcat('innerSolarSystem'),'png');
    elseif(j==5)
        saveas(gcf,strcat('barycentricInnerSolarSystem'),'png');
    elseif(j==6)
        saveas(gcf,strcat('hierarchicalInnerSolarSystem'),'png');
    elseif(j==7)
        saveas(gcf,strcat('perturbationInfluence'),'png');
    elseif(j==8)
        saveas(gcf,strcat('thrustAlongVelocityVectorResults'),'png');
    elseif(j==9)
        saveas(gcf,strcat('thrustFromFileResults'),'png');
    elseif(j==10)
        saveas(gcf,strcat('thrustFromFileProfiles'),'png');
    end

end
