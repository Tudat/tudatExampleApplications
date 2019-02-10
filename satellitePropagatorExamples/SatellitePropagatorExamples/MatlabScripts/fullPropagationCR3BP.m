clear all;
clc;
close all;

propagatedStateHistory = importdata('../SimulationOutput/fullProblemPropagation.dat');
normalisedPropagatedStateHistory = importdata('../SimulationOutput/fullProblemPropagationNormalisedCoRotatingFrame.dat');
CR3BPstateHistory = importdata('../SimulationOutput/CR3BPsolution.dat');
CR3BPunnormalisedStateHistory = importdata('../SimulationOutput/CR3BPnormalisedCoRotatingFrame.dat');

propagatedStateHistoryPerturbedCase = importdata('../SimulationOutput/fullProblemPropagationPerturbedCase.dat');
normalisedPropagatedStateHistoryPerturbedCase = importdata('../SimulationOutput/fullProblemPropagationNormalisedCoRotatingFramePerturbedCase.dat');
CR3BPstateHistoryPerturbedCase = importdata('../SimulationOutput/CR3BPsolutionPerturbedCase.dat');
CR3BPunnormalisedStateHistoryPerturbedCase = importdata('../SimulationOutput/CR3BPnormalisedCoRotatingFramePerturbedCase.dat');

% Difference in cartesian position coordinates between CR3BP and propagated state 
%  (unperturbed case respecting the CR3BP assumptions)
difference_x = CR3BPstateHistory(1:end,2) - propagatedStateHistory(:,2);
difference_y = CR3BPstateHistory(1:end,3) - propagatedStateHistory(:,3);
difference_z = CR3BPstateHistory(1:end,4) - propagatedStateHistory(:,4);

% Difference in cartesian position coordinates between CR3BP and propagated state (perturbed case) 
difference_x_perturbed_case = CR3BPstateHistoryPerturbedCase(1:end,2) - propagatedStateHistoryPerturbedCase(:,2);
difference_y_perturbed_case = CR3BPstateHistoryPerturbedCase(1:end,3) - propagatedStateHistoryPerturbedCase(:,3);
difference_z_perturbed_case = CR3BPstateHistoryPerturbedCase(1:end,4) - propagatedStateHistoryPerturbedCase(:,4);

% Difference in state position for the unperturbed case
figure();
hold on;
plot(propagatedStateHistory(:,1), difference_x);
plot(propagatedStateHistory(:,1), difference_y);
plot(propagatedStateHistory(:,1), difference_z);
hold off;
grid;
ylabel('Difference in position coordinates [m]');
xlabel('Time [s]');
legend({'x', 'y', 'z'});

% Difference in state position for the perturbed case
figure();
hold on;
plot(propagatedStateHistoryPerturbedCase(:,1), difference_x_perturbed_case);
plot(propagatedStateHistoryPerturbedCase(:,1), difference_y_perturbed_case);
plot(propagatedStateHistoryPerturbedCase(:,1), difference_z_perturbed_case);
hold off;
grid;
ylabel('Difference in position coordinates [m]');
xlabel('Time [s]');
legend({'x', 'y', 'z'});

% Spacecraft trajectory in the rotating normalised reference frame (unperturbed case)
figure();
hold on;
plot(normalisedPropagatedStateHistory(:,2), normalisedPropagatedStateHistory(:,3));
plot(CR3BPunnormalisedStateHistory(:,2), CR3BPunnormalisedStateHistory(:,3));
hold off;
grid;
xlabel('x [m]');
ylabel('y [m]');
xlim([-1.5 1.5]);
legend({'Full problem propagation', 'CR3BP solution'});

% Spacecraft trajectory in the rotating normalised reference frame (perturbed case)
figure();
hold on;
plot(normalisedPropagatedStateHistoryPerturbedCase(:,2), normalisedPropagatedStateHistoryPerturbedCase(:,3));
plot(CR3BPunnormalisedStateHistoryPerturbedCase(:,2), CR3BPunnormalisedStateHistoryPerturbedCase(:,3));
hold off;
grid;
xlabel('x [m]');
ylabel('y [m]');
legend({'Full problem propagation', 'CR3BP solution'});
