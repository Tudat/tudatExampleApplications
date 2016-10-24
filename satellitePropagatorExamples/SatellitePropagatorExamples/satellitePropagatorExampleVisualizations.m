clc
close all
clear all
folder = '';
%%

figure(1)

load(strcat(folder,'asterixPropagationHistory.dat'));
load(strcat(folder,'obelixPropagationHistory.dat'));
scatter3(asterixPropagationHistory(:,2),asterixPropagationHistory(:,3),asterixPropagationHistory(:,4))
hold on
scatter3(obelixPropagationHistory(:,2),obelixPropagationHistory(:,3),obelixPropagationHistory(:,4))
axis equal

legend('Asterix orbit','Obelix orbit')

xlabel('x_{J2000} [m]')
ylabel('y_{J2000} [m]')
zlabel('z_{J2000} [m]')
%%
figure(2)

for i=1:30
    currentGalileo = load(strcat(folder,'galileoSatellite',num2str(i),'.dat'));
    scatter3(currentGalileo(:,2),currentGalileo(:,3),currentGalileo(:,4))
    hold on
    axis equal
end
xlabel('x_{J2000} [m]')
ylabel('y_{J2000} [m]')
zlabel('z_{J2000} [m]')

%%
figure(3)

load(strcat(folder,'apolloPropagationHistory.dat'));

subplot(1,2,1)
distance = (sqrt(apolloPropagationHistory(:,2).^2+apolloPropagationHistory(:,3).^2+apolloPropagationHistory(:,4).^2) - 6378E3)/1000;
plot(apolloPropagationHistory(:,1),distance);
grid on
xlabel('Time since propagation start [s]')
ylabel('Earth altitude [km]')

subplot(1,2,2)
speed = sqrt(apolloPropagationHistory(:,5).^2+apolloPropagationHistory(:,6).^2+apolloPropagationHistory(:,7).^2);
plot(apolloPropagationHistory(:,1),speed);
grid on
xlabel('Time since propagation start [s]')
ylabel('Earth-centered inertial speed [m/s]')

%%

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
    barycentricSolarSystemOrbits{i} = load(strcat('innerSolarSystemPropagationHistory',propagatedBodies{i},'0.dat'));
    figure(4)
    scatter3(barycentricSolarSystemOrbits{i}(:,2),barycentricSolarSystemOrbits{i}(:,3),barycentricSolarSystemOrbits{i}(:,4))
    hold on
    
    figure(5)
    subplot(3,2,i)
    scatter3(barycentricSolarSystemOrbits{i}(:,2),barycentricSolarSystemOrbits{i}(:,3),barycentricSolarSystemOrbits{i}(:,4))
    title(strcat('Barycentric orbit of ',propagatedBodies{i}));
    
    localBarycentricSolarSystemOrbits{i} = load(strcat('innerSolarSystemPropagationHistory',propagatedBodies{i},'1.dat'));
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

figure(7)
perturbedSatellite = load('singlePerturbedSatellitePropagationHistory.dat');
unperturbedSatellite = load('singleSatellitePropagationHistory.dat');

for i=1:3
subplot(3,2,2*(i-1)+1)
plot(perturbedSatellite(:,1),perturbedSatellite(:,i+1))
hold on
plot(unperturbedSatellite(:,1),unperturbedSatellite(:,i+1),'r--')
subplot(3,2,2*i)
plot(perturbedSatellite(:,1),perturbedSatellite(:,i+1)-unperturbedSatellite(:,i+1))

end

%%
close all

load('thrustPropagationHistory.dat')
plot3(thrustPropagationHistory(:,2),thrustPropagationHistory(:,3),thrustPropagationHistory(:,4))
axis equal
