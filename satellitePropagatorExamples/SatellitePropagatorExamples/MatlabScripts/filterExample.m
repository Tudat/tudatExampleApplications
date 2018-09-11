fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;

%...Main path
directory = '../SimulationOutput/FilterEstimation/';

%% Load Actual State and Noise

%...Actual state
filename = fullfile(directory,'actualStateHistory.dat');
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f %f','CollectOutput',true,'Delimiter',',');
actualStateHistory = cpp_result{1}(:,2:end)';
fclose(fileID);

%...System noise
filename = fullfile(directory,'systemNoise.dat');
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f','CollectOutput',true);
systemNoiseHistory = cpp_result{1}';
fclose(fileID);

%...Measurement noise
filename = fullfile(directory,'measurementNoise.dat');
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f','CollectOutput',true);
measurementNoiseHistory = cpp_result{1}';
fclose(fileID);

%% Load C++ EKF Output

%...EKF state
filename = fullfile(directory,'EKFEstimatedStateHistory.dat');
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f %f','CollectOutput',true,'Delimiter',',');
ekfStateHistory = cpp_result{1}(:,2:end)';
fclose(fileID);

%...EKF covariance
filename = fullfile(directory,'EKFEstimatedCovarianceHistory.dat');
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f','CollectOutput',true,'Delimiter',',');
ekfCovarianceHistory = cpp_result{1}(:,2:end)';
fclose(fileID);

%% Load C++ UKF Output

%...UKF state
filename = fullfile(directory,'UKFEstimatedStateHistory.dat');
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f %f','CollectOutput',true,'Delimiter',',');
ukfStateHistory = cpp_result{1}(:,2:end)';
fclose(fileID);

%...UKF covariance
filename = fullfile(directory,'UKFEstimatedCovarianceHistory.dat');
fileID = fopen(filename,'r');
cpp_result = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f','CollectOutput',true,'Delimiter',',');
ukfCovarianceHistory = cpp_result{1}(:,2:end)';
fclose(fileID);

%% Plot States

xAxisValue = 1:size(actualStateHistory,2);

%...Plot states over time
figure

%...Position
subplot(1,3,1)
hold on
plot(xAxisValue,ekfStateHistory(1,:),'LineWidth',1.25)
plot(xAxisValue,ukfStateHistory(1,:),'LineWidth',1.25)
plot(xAxisValue,actualStateHistory(1,:),'LineWidth',1.25)
hold off
grid on
legend('EKF','UKF','Actual','Location','Best')
set(gca,'FontSize',15)
title('x_1')

%...Velocity
subplot(1,3,2)
hold on
plot(xAxisValue,ekfStateHistory(2,:),'LineWidth',1.25)
plot(xAxisValue,ukfStateHistory(2,:),'LineWidth',1.25)
plot(xAxisValue,actualStateHistory(2,:),'LineWidth',1.25)
hold off
grid on
legend('EKF','UKF','Actual','Location','Best')
set(gca,'FontSize',15)
title('x_2')

%...Ballistic coefficient
subplot(1,3,3)
hold on
plot(xAxisValue,ekfStateHistory(3,:),'LineWidth',1.25)
plot(xAxisValue,ukfStateHistory(3,:),'LineWidth',1.25)
plot(xAxisValue,actualStateHistory(3,:),'LineWidth',1.25)
hold off
grid on
legend('EKF','UKF','Actual','Location','Best')
set(gca,'FontSize',15)
title('x_3')

%% Plot State Differences

%...Plot state errors over time
F = figure( 'rend', 'painters', 'pos', [ 440, 378, 840, 630 ] );

%...Position
subplot(1,3,1)
hold on
plot(xAxisValue,ekfStateHistory(1,:)-actualStateHistory(1,:),'LineWidth',1.25)
plot(xAxisValue,ukfStateHistory(1,:)-actualStateHistory(1,:),'LineWidth',1.25)
plot(xAxisValue,sqrt(ekfCovarianceHistory(1,:)),'LineWidth',1.25,'LineStyle','--','Color',[0.929,0.694,0.125])
plot(xAxisValue,sqrt(ukfCovarianceHistory(1,:)),'LineWidth',1.25,'LineStyle','-.','Color',[0.494,0.184,0.556])
plot(xAxisValue,-sqrt(ekfCovarianceHistory(1,:)),'LineWidth',1.25,'LineStyle','--','Color',[0.929,0.694,0.125])
plot(xAxisValue,-sqrt(ukfCovarianceHistory(1,:)),'LineWidth',1.25,'LineStyle','-.','Color',[0.494,0.184,0.556])
hold off
grid on
legend('EKF','UKF','STD EKF','STD UKF','Location','Best')
set(gca,'FontSize',15)
title('x_1')

%...Velocity
subplot(1,3,2)
hold on
plot(xAxisValue,ekfStateHistory(2,:)-actualStateHistory(2,:),'LineWidth',1.25)
plot(xAxisValue,ukfStateHistory(2,:)-actualStateHistory(2,:),'LineWidth',1.25)
plot(xAxisValue,sqrt(ekfCovarianceHistory(5,:)),'LineWidth',1.25,'LineStyle','--','Color',[0.929,0.694,0.125])
plot(xAxisValue,sqrt(ukfCovarianceHistory(5,:)),'LineWidth',1.25,'LineStyle','-.','Color',[0.494,0.184,0.556])
plot(xAxisValue,-sqrt(ekfCovarianceHistory(5,:)),'LineWidth',1.25,'LineStyle','--','Color',[0.929,0.694,0.125])
plot(xAxisValue,-sqrt(ukfCovarianceHistory(5,:)),'LineWidth',1.25,'LineStyle','-.','Color',[0.494,0.184,0.556])
hold off
grid on
legend('EKF','UKF','STD EKF','STD UKF','Location','Best')
set(gca,'FontSize',15)
title('x_2')

%...Ballistic coefficient
subplot(1,3,3)
hold on
plot(xAxisValue,ekfStateHistory(3,:)-actualStateHistory(3,:),'LineWidth',1.25)
plot(xAxisValue,ukfStateHistory(3,:)-actualStateHistory(3,:),'LineWidth',1.25)
plot(xAxisValue,sqrt(ekfCovarianceHistory(9,:)),'LineWidth',1.25,'LineStyle','--','Color',[0.929,0.694,0.125])
plot(xAxisValue,sqrt(ukfCovarianceHistory(9,:)),'LineWidth',1.25,'LineStyle','-.','Color',[0.494,0.184,0.556])
plot(xAxisValue,-sqrt(ekfCovarianceHistory(9,:)),'LineWidth',1.25,'LineStyle','--','Color',[0.929,0.694,0.125])
plot(xAxisValue,-sqrt(ukfCovarianceHistory(9,:)),'LineWidth',1.25,'LineStyle','-.','Color',[0.494,0.184,0.556])
hold off
grid on
legend('EKF','UKF','STD EKF','STD UKF','Location','Best')
set(gca,'FontSize',15)
title('x_3')

%...Save figure
saveas( F, 'filterExample', 'epsc' )