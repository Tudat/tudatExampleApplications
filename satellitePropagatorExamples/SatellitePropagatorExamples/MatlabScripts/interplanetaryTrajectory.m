clc
close all
clear all

%%
%%%%% DEFINE DATA DIRECTORY AND LOAD FILES
dataDirectory = '../SimulationOutput/';

earthTrajectory = load(strcat(dataDirectory,'earthTrajectory.dat'));
venusTrajectory = load(strcat(dataDirectory,'venusTrajectory.dat'));
mercuryTrajectory = load(strcat(dataDirectory,'mercuryTrajectory.dat'));

messengerTrajectory = load(strcat(dataDirectory,'messengerTrajectory.dat'));
messengerManeuvers = load(strcat(dataDirectory,'messengerManeuvers.dat'));

%%

figure;
plot(messengerTrajectory(:,2),messengerTrajectory(:,3));
hold on
plot(earthTrajectory(:,2),earthTrajectory(:,3),'--');
plot(venusTrajectory(:,2),venusTrajectory(:,3),'--');
plot(mercuryTrajectory(:,2),mercuryTrajectory(:,3),'--');
scatter(messengerManeuvers(:,2),messengerManeuvers(:,3), 'Filled');
grid on
xlabel('x[m]');
ylabel('y[m]');
zlabel('z[m]');