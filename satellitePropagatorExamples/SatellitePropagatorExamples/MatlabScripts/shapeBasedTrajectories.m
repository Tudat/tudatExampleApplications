%
% This script processes the results of the shape-based trajectories example, 
% run by the shapeBasedTrajectoryDesign.cpp satellitePropagatorExamples application.
%

set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

clc
close all
clear all

%% Load results

dataDirectory = '../SimulationOutput/ShapeBasedTrajectoriesExample/';

% Retrieve shape-based trajectories
hodographicallyShapedTrajectory = load(strcat(dataDirectory,'hodographicShapingAnalyticalResults.dat'));
sphericallyShapedTrajectory = load(strcat(dataDirectory,'sphericalShapingAnalyticalResults.dat'));

% Retrieve results of the full propagation for the unperturbed case
hodographicShapingPropagationUnperturbedCase = load(strcat(dataDirectory,'hodographicShapingPropagationUnperturbedCase.dat'));
sphericalShapingPropagationUnperturbedCase = load(strcat(dataDirectory,'sphericalShapingPropagationUnperturbedCase.dat'));

% Retrieve results of the full propagation for the perturbed case
hodographicShapingPropagationPerturbedCase = load(strcat(dataDirectory,'hodographicShapingPropagationPerturbedCase.dat'));
sphericalShapingPropagationPerturbedCase = load(strcat(dataDirectory,'sphericalShapingPropagationPerturbedCase.dat'));

% Retrieve mass, thrust, and thrust acceleration profiles for hodographic shaping
hodographicShapingMassProfile = load(strcat(dataDirectory,'hodographicShapingMassProfile.dat'));
hodographicShapingThrustProfile = load(strcat(dataDirectory,'hodographicShapingThrustProfile.dat'));
hodographicShapingThrustAccelerationProfile = load(strcat(dataDirectory,'hodographicShapingThrustAccelerationProfile.dat'));

% Retrieve mass, thrust, and thrust acceleration profiles for spherical shaping
sphericalShapingMassProfile = load(strcat(dataDirectory,'sphericalShapingMassProfile.dat'));
sphericalShapingThrustProfile = load(strcat(dataDirectory,'sphericalShapingThrustProfile.dat'));
sphericalShapingThrustAccelerationProfile = load(strcat(dataDirectory,'sphericalShapingThrustAccelerationProfile.dat'));


%% Plots analytical results from shape-based methods

figure(1);

subplot(2,2,1);
% Plot the hodographic and spherical shaping trajectories
hold on;
plot3(hodographicallyShapedTrajectory(:,2),hodographicallyShapedTrajectory(:,3), ...
    hodographicallyShapedTrajectory(:,4), 'r', 'lineWidth', 1.8);
plot3(sphericallyShapedTrajectory(:,2),sphericallyShapedTrajectory(:,3), ...
    sphericallyShapedTrajectory(:,4), 'b', 'lineWidth', 1.8);
hold off;
title('Spacecraft trajectory', 'fontSize', 12);
legend({'Hodographic shaping', 'Spherical shaping'}, 'fontSize', 12, 'location', 'Best');
grid;
xlabel('x[m]', 'fontSize', 12);
ylabel('y[m]', 'fontSize', 12);
zlabel('z[m]', 'fontSize', 12);
view([30,15]);

subplot(2,2,2);
% Plot the thrust acceleration profiles for both hodographic and spherical shaping
hold on;
plot(hodographicShapingThrustAccelerationProfile(:,1),hodographicShapingThrustAccelerationProfile(:,2), ...
    'r', 'lineWidth', 1.8);
plot(sphericalShapingThrustAccelerationProfile(:,1),sphericalShapingThrustAccelerationProfile(:,2), ...
    'b', 'lineWidth', 1.8);
hold off;
title('Acceleration profile', 'fontSize', 12);
grid;
xlabel('Time [s]', 'fontSize', 12);
ylabel('Thrust acceleration [$\mathrm{ms^{-2}}$]', 'fontSize', 12);
xlim([0 hodographicallyShapedTrajectory(end,1)]);

subplot(2,2,3);
% Plot the mass profiles for both hodographic and spherical shaping
hold on;
plot(hodographicShapingMassProfile(:,1),hodographicShapingMassProfile(:,2), 'r', 'lineWidth', 1.8);
plot( sphericalShapingMassProfile(:,1),sphericalShapingMassProfile(:,2), 'b', 'lineWidth', 1.8);
hold off;
title('Mass profile', 'fontSize', 12);
grid;
xlabel('Time [s]', 'fontSize', 12);
ylabel('Mass [kg]', 'fontSize', 12);
xlim([0 hodographicallyShapedTrajectory(end,1)]);

subplot(2,2,4);
% Plot the thrust profiles for both hodographic and spherical shaping
hold on;
plot(hodographicShapingThrustProfile(:,1),hodographicShapingThrustProfile(:,2), 'r', 'lineWidth', 1.8);
plot( sphericalShapingThrustProfile(:,1),sphericalShapingThrustProfile(:,2), 'b', 'lineWidth', 1.8);
hold off;
title('Thrust profile', 'fontSize', 12);
grid;
xlabel('Time [s]', 'fontSize', 12);
ylabel('Thrust [N]', 'fontSize', 12);
xlim([0 hodographicallyShapedTrajectory(end,1)]);


%% Plots results of the full propagation (unperturbed and perturbed cases)

figure(2);

subplot(2,2,1);
% Plot difference between analytical and full propagation results for
% hodographic shaping (unperturbed case)
plot(hodographicallyShapedTrajectory(:,1), ...
    hodographicallyShapedTrajectory(:,2) - hodographicShapingPropagationUnperturbedCase(:,2), 'lineWidth', 1.8 );
hold on
plot(hodographicallyShapedTrajectory(:,1), ...
    hodographicallyShapedTrajectory(:,3) - hodographicShapingPropagationUnperturbedCase(:,3), 'lineWidth', 1.8 );
plot(hodographicallyShapedTrajectory(:,1), ...
    hodographicallyShapedTrajectory(:,4) - hodographicShapingPropagationUnperturbedCase(:,4), 'lineWidth', 1.8 );
ylabel('Difference analytical/numerical [m]', 'fontSize', 12);
xlabel('Time [s]', 'fontSize', 12);
grid;
legend({'x [m]', 'y [m]', 'z [m]'}, 'fontSize', 10, 'location', 'Best', 'fontSize', 12);
xlim([0 hodographicallyShapedTrajectory(end,1)]);
title('Hodographic shaping - simplified case', 'fontSize', 12);

subplot(2,2,2);
% Plot difference between analytical and full propagation results for
% spherical shaping (unperturbed case)
plot(sphericallyShapedTrajectory(10:end-10,1), ...
    sphericallyShapedTrajectory(10:end-10,2) - sphericalShapingPropagationUnperturbedCase(10:end-10,2), ...
    'lineWidth', 1.8 );
hold on
plot(sphericallyShapedTrajectory(10:end-10,1), ...
    sphericallyShapedTrajectory(10:end-10,3) - sphericalShapingPropagationUnperturbedCase(10:end-10,3), ...
    'lineWidth', 1.8 );
plot(sphericallyShapedTrajectory(10:end-10,1), ...
    sphericallyShapedTrajectory(10:end-10,4) - sphericalShapingPropagationUnperturbedCase(10:end-10,4), ...
    'lineWidth', 1.8 );
ylabel('Difference analytical/numerical [m]', 'fontSize', 12);
xlabel('Time [s]', 'fontSize', 12);
grid;
legend({'x [m]', 'y [m]', 'z [m]'}, 'fontSize', 12, 'location', 'Best');
xlim([sphericallyShapedTrajectory(10,1) sphericallyShapedTrajectory(end-10,1)]);
title('Spherical shaping - simplified case', 'fontSize', 12);

subplot(2,2,3);
% Plot difference between analytical and full propagation results for
% hodographic shaping (perturbed case)
plot(hodographicallyShapedTrajectory(:,1), ...
    hodographicallyShapedTrajectory(:,2) - hodographicShapingPropagationPerturbedCase(:,2), 'lineWidth', 1.8 );
hold on
plot(hodographicallyShapedTrajectory(:,1), ...
    hodographicallyShapedTrajectory(:,3) - hodographicShapingPropagationPerturbedCase(:,3), 'lineWidth', 1.8 );
plot(hodographicallyShapedTrajectory(:,1), ...
    hodographicallyShapedTrajectory(:,4) - hodographicShapingPropagationPerturbedCase(:,4), 'lineWidth', 1.8 );
ylabel('Difference analytical/numerical [m]', 'fontSize', 12);
xlabel('Time [s]', 'fontSize', 12);
grid;
legend({'x [m]', 'y [m]', 'z [m]'}, 'fontSize', 12, 'location', 'Best');
xlim([0 hodographicallyShapedTrajectory(end,1)]);
title('Hodographic shaping - perturbed case', 'fontSize', 12);

subplot(2,2,4);
% Plot difference between analytical and full propagation results for
% spherical shaping (perturbed case)
plot(sphericallyShapedTrajectory(:,1), ...
    sphericallyShapedTrajectory(:,2) - sphericalShapingPropagationPerturbedCase(:,2), 'lineWidth', 1.8 );
hold on
plot(sphericallyShapedTrajectory(:,1), ...
    sphericallyShapedTrajectory(:,3) - sphericalShapingPropagationPerturbedCase(:,3), 'lineWidth', 1.8 );
plot(sphericallyShapedTrajectory(:,1), ...
    sphericallyShapedTrajectory(:,4) - sphericalShapingPropagationPerturbedCase(:,4), 'lineWidth', 1.8 );
ylabel('Difference analytical/numerical [m]', 'fontSize', 12);
xlabel('Time [s]', 'fontSize', 12);
grid;
legend({'x [m]', 'y [m]', 'z [m]'}, 'fontSize', 12, 'location', 'Best');
xlim([0 sphericallyShapedTrajectory(end,1)]);
title('Spherical shaping - perturbed case', 'fontSize', 12);


