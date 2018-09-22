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

xAxisValue = cpp_result{1}(:,1);

%...Plot variables
labels = {'Position [ft]','Velocity [ft s^{-1}]','Ballistic Coefficient [lb ft^{-2}]'};
titles = {'x_1','x_2','x_3'};

%...Plot states over time
figure( 'rend', 'painters', 'pos', [ 440, 378, 840, 630 ] )
for i = 1:3
    subplot(1,3,i)
    hold on
    plot(xAxisValue,ekfStateHistory(i,:),'LineWidth',1.25)
    plot(xAxisValue,ukfStateHistory(i,:),'LineWidth',1.25,'LineStyle','--')
    plot(xAxisValue,actualStateHistory(i,:),'LineWidth',1.25,'LineStyle','-.')
    hold off
    grid on
    xlabel('Time [s]')
    ylabel(labels{i})
    legend('EKF','UKF','Actual','Location','Best')
    set(gca,'FontSize',15)
    title(titles{i})
end

%% Plot State Errors

%...Plot variables
labels = {'Position Error [ft]','Velocity Error [ft s^{-1}]','Ballistic Coefficient Error [lb ft^{-2}]'};
titles = {'x_1','x_2','x_3'};

%...Plot state errors over time
F = figure( 'rend', 'painters', 'pos', [ 440, 378, 840, 630 ] );
for i = 1:3
    subplot(1,3,i)
    hold on
    plot(xAxisValue,ekfStateHistory(i,:)-actualStateHistory(i,:),'LineWidth',1.25)
    plot(xAxisValue,ukfStateHistory(i,:)-actualStateHistory(i,:),'LineWidth',1.25,'LineStyle','--')
    if i == 1, j = 1; elseif i == 2, j = 5; else, j = 9; end
    plot(xAxisValue,sqrt(ekfCovarianceHistory(j,:)),'LineWidth',1.25,'LineStyle','-.','Color',[0.929,0.694,0.125])
    plot(xAxisValue,sqrt(ukfCovarianceHistory(j,:)),'LineWidth',1.25,'LineStyle',':','Color',[0.494,0.184,0.556])
    plot(xAxisValue,-sqrt(ekfCovarianceHistory(j,:)),'LineWidth',1.25,'LineStyle','-.','Color',[0.929,0.694,0.125])
    plot(xAxisValue,-sqrt(ukfCovarianceHistory(j,:)),'LineWidth',1.25,'LineStyle',':','Color',[0.494,0.184,0.556])
    hold off
    grid on
    xlabel('Time [s]')
    ylabel(labels{i})
    legend('EKF','UKF','STD EKF','STD UKF','Location','Best')
    set(gca,'FontSize',15)
    title(titles{i})
end

%...Save figure
saveas( F, 'filterExample', 'png' )