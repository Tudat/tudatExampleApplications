
clc
close all
clear all

% Load hodographic shaping results (trajectory, mass, thrust, and thrust acceleration profiles)
hodographicShapingTrajectory = load('../SimulationOutput/hodographicShapingTrajectory.dat');
hodographicShapingMassProfile = load('../SimulationOutput/hodographicShapingMassProfile.dat');
hodographicShapingThrustProfile = load('../SimulationOutput/hodographicShapingThrustProfile.dat');
hodographicShapingThrustAccelerationProfile = load('../SimulationOutput/hodographicShapingThrustAcceleration.dat');

figure

plot3( hodographicShapingTrajectory(:,2), hodographicShapingTrajectory(:,3), hodographicShapingTrajectory(:,4), 'b', 'lineWidth', 1.5 );
grid on

figure

for i=1:3
    subplot(1,3,1)
    plot(hodographicShapingThrustAccelerationProfile(:,1),hodographicShapingThrustAccelerationProfile(:,i+1))
    hold on
    grid on
    
end

subplot(1,3,2)
plot(sqrt(sum(hodographicShapingThrustProfile(:,2:4).^2')))
grid on

subplot(1,3,3)
plot(hodographicShapingMassProfile(:,1),hodographicShapingMassProfile(:,2))
grid on

%%
for i=1:5
    hodographicMultiObjectiveFitness = load(strcat('../SimulationOutput/fitness_hodographic_multi',num2str(i),'.dat'));
    
    
    scatter(hodographicMultiObjectiveFitness(:,1)/1000,hodographicMultiObjectiveFitness(:,2))
    hold on
end
grid on
xlabel('Delta V [m/s]')
ylabel('Maximum acceleration [m/s2]')
