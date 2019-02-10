clear all; 
clc;
close all;

dataDirectory = '../SimulationOutput/';
earthTrajectory = load(strcat(dataDirectory,'earthTrajectory.dat'));
venusTrajectory = load(strcat(dataDirectory,'venusTrajectory.dat'));
mercuryTrajectory = load(strcat(dataDirectory,'mercuryTrajectory.dat'));

% Ideal case: propagation of the spacecraft under the patched conics assumptions

propagatedStateDepartureBeforeDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_0_leg_0.dat');
propagatedStateDepartureAfterDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_0_leg_1.dat');
propagatedStateSwingby1BeforeDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_0_leg_2.dat');
propagatedStateSwingby1AfterDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_0_leg_3.dat');
propagatedStateSwingby2BeforeDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_0_leg_4.dat');
propagatedStateSwingby2AfterDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_0_leg_5.dat');
propagatedStateSwingby3BeforeDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_0_leg_6.dat');
propagatedStateSwingby3AfterDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_0_leg_7.dat');

patchedConicsDepartureBeforeDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_0_leg_0.dat');
patchedConicsDepartureAfterDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_0_leg_1.dat');
patchedConicsSwingby1BeforeDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_0_leg_2.dat');
patchedConicsSwingby1AfterDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_0_leg_3.dat');
patchedConicsSwingby2BeforeDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_0_leg_4.dat');
patchedConicsSwingby2AfterDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_0_leg_5.dat');
patchedConicsSwingby3BeforeDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_0_leg_6.dat');
patchedConicsSwingby3AfterDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_0_leg_7.dat');

differenceDepartureLegBeforeDsm = [propagatedStateDepartureBeforeDsm(:,1) sqrt( (propagatedStateDepartureBeforeDsm(:,2) - patchedConicsDepartureBeforeDsm(:,2)).^2 + (propagatedStateDepartureBeforeDsm(:,3) - patchedConicsDepartureBeforeDsm(:,3)).^2 + (propagatedStateDepartureBeforeDsm(:,4) - patchedConicsDepartureBeforeDsm(:,4)).^2)];
differenceDepartureLegAfterDsm = [propagatedStateDepartureAfterDsm(:,1) sqrt( (propagatedStateDepartureAfterDsm(:,2) - patchedConicsDepartureAfterDsm(:,2)).^2 + (propagatedStateDepartureAfterDsm(:,3) - patchedConicsDepartureAfterDsm(:,3)).^2 + (propagatedStateDepartureAfterDsm(:,4) - patchedConicsDepartureAfterDsm(:,4)).^2)];
differenceSwingby1BeforeDsm = [propagatedStateSwingby1BeforeDsm(:,1) sqrt( (propagatedStateSwingby1BeforeDsm(:,2) - patchedConicsSwingby1BeforeDsm(:,2)).^2 + (propagatedStateSwingby1BeforeDsm(:,3) - patchedConicsSwingby1BeforeDsm(:,3)).^2 + (propagatedStateSwingby1BeforeDsm(:,4) - patchedConicsSwingby1BeforeDsm(:,4)).^2)];
differenceSwingby1AfterDsm = [propagatedStateSwingby1AfterDsm(:,1) sqrt( (propagatedStateSwingby1AfterDsm(:,2) - patchedConicsSwingby1AfterDsm(:,2)).^2 + (propagatedStateSwingby1AfterDsm(:,3) - patchedConicsSwingby1AfterDsm(:,3)).^2 + (propagatedStateSwingby1AfterDsm(:,4) - patchedConicsSwingby1AfterDsm(:,4)).^2)];
differenceSwingby2BeforeDsm = [propagatedStateSwingby2BeforeDsm(:,1) sqrt( (propagatedStateSwingby2BeforeDsm(:,2) - patchedConicsSwingby2BeforeDsm(:,2)).^2 + (propagatedStateSwingby2BeforeDsm(:,3) - patchedConicsSwingby2BeforeDsm(:,3)).^2 + (propagatedStateSwingby2BeforeDsm(:,4) - patchedConicsSwingby2BeforeDsm(:,4)).^2)];
differenceSwingby2AfterDsm = [propagatedStateSwingby2AfterDsm(:,1) sqrt( (propagatedStateSwingby2AfterDsm(:,2) - patchedConicsSwingby2AfterDsm(:,2)).^2 + (propagatedStateSwingby2AfterDsm(:,3) - patchedConicsSwingby2AfterDsm(:,3)).^2 + (propagatedStateSwingby2AfterDsm(:,4) - patchedConicsSwingby2AfterDsm(:,4)).^2)];
differenceSwingby3BeforeDsm = [propagatedStateSwingby3BeforeDsm(:,1) sqrt( (propagatedStateSwingby3BeforeDsm(:,2) - patchedConicsSwingby3BeforeDsm(:,2)).^2 + (propagatedStateSwingby3BeforeDsm(:,3) - patchedConicsSwingby3BeforeDsm(:,3)).^2 + (propagatedStateSwingby3BeforeDsm(:,4) - patchedConicsSwingby3BeforeDsm(:,4)).^2)];
differenceSwingby3AfterDsm = [propagatedStateSwingby3AfterDsm(:,1) sqrt( (propagatedStateSwingby3AfterDsm(:,2) - patchedConicsSwingby3AfterDsm(:,2)).^2 + (propagatedStateSwingby3AfterDsm(:,3) - patchedConicsSwingby3AfterDsm(:,3)).^2 + (propagatedStateSwingby3AfterDsm(:,4) - patchedConicsSwingby3AfterDsm(:,4)).^2)];


figure();
hold on;
% Orbits of the Earth, Venus and Mercury
plot(earthTrajectory(:,2), earthTrajectory(:,3), '--k');
plot(venusTrajectory(:,2), venusTrajectory(:,3), '--k');
plot(mercuryTrajectory(:,2), mercuryTrajectory(:,3), '--k');

% Propagation of the spacecraft along the Messenger trajectory (under the patched conics assumptions)
plot(propagatedStateDepartureBeforeDsm(:,2), propagatedStateDepartureBeforeDsm(:,3), 'color', [0, 0.4470, 0.7410], 'lineWidth', 1.1);
plot(propagatedStateDepartureAfterDsm(:,2), propagatedStateDepartureAfterDsm(:,3), 'color', [0, 0.4470, 0.7410], 'lineWidth', 1.1);
plot(propagatedStateSwingby1BeforeDsm(:,2), propagatedStateSwingby1BeforeDsm(:,3),'color', [0.8500, 0.3250, 0.0980], 'lineWidth', 1.1);
plot(propagatedStateSwingby1AfterDsm(:,2), propagatedStateSwingby1AfterDsm(:,3), 'color', [0.8500, 0.3250, 0.0980], 'lineWidth', 1.1);
plot(propagatedStateSwingby2BeforeDsm(:,2), propagatedStateSwingby2BeforeDsm(:,3), 'color', [0.9290, 0.6940, 0.1250], 'lineWidth', 1.1);
plot(propagatedStateSwingby2AfterDsm(:,2), propagatedStateSwingby2AfterDsm(:,3), 'color', [0.9290, 0.6940, 0.1250], 'lineWidth', 1.1);
plot(propagatedStateSwingby3BeforeDsm(:,2), propagatedStateSwingby3BeforeDsm(:,3), 'color', [0.4940, 0.1840, 0.5560], 'lineWidth', 1.1);
plot(propagatedStateSwingby3AfterDsm(:,2), propagatedStateSwingby3AfterDsm(:,3), 'color', [0.4940, 0.1840, 0.5560], 'lineWidth', 1.1);

% Transfer bodies of the trajectory
plot(patchedConicsDepartureBeforeDsm(1,2), patchedConicsDepartureBeforeDsm(1,3), '.k', 'markerSize', 12);
plot(patchedConicsSwingby1BeforeDsm(1,2), patchedConicsSwingby1BeforeDsm(1,3), '.k', 'markerSize', 12);
plot(patchedConicsSwingby2BeforeDsm(1,2), patchedConicsSwingby2BeforeDsm(1,3), '.k', 'markerSize', 12);
plot(patchedConicsSwingby3BeforeDsm(1,2), patchedConicsSwingby3BeforeDsm(1,3), '.k', 'markerSize', 12);
plot(patchedConicsSwingby3AfterDsm(end,2), patchedConicsSwingby3AfterDsm(end,3), '.k', 'markerSize', 12)

% DSM locations along the trajectory
plot(patchedConicsDepartureBeforeDsm(end,2), patchedConicsDepartureBeforeDsm(end,3), '.r', 'markerSize', 12);
plot(patchedConicsSwingby1BeforeDsm(end,2), patchedConicsSwingby1BeforeDsm(end,3), '.r', 'markerSize', 12);
plot(patchedConicsSwingby2BeforeDsm(end,2), patchedConicsSwingby2BeforeDsm(end,3), '.r', 'markerSize', 12);
plot(patchedConicsSwingby3BeforeDsm(end,2), patchedConicsSwingby3BeforeDsm(end,3), '.r', 'markerSize', 12);

hold off;
grid;
xlabel('x [m]', 'fontSize', 12 );
ylabel('y [m]', 'fontSize', 12 );
title('Full problem under the patched conics assumptions');



figure();
semilogy(differenceDepartureLegBeforeDsm(:,1), differenceDepartureLegBeforeDsm(:,2), 'color', [0, 0.4470, 0.7410]);
hold on;
plot(differenceDepartureLegAfterDsm(:,1), differenceDepartureLegAfterDsm(:,2),'color', [0, 0.4470, 0.7410]);
plot(differenceSwingby1BeforeDsm(:,1), differenceSwingby1BeforeDsm(:,2),'color', [0.8500, 0.3250, 0.0980]);
plot(differenceSwingby1AfterDsm(:,1), differenceSwingby1AfterDsm(:,2),'color', [0.8500, 0.3250, 0.0980]);
plot(differenceSwingby2BeforeDsm(:,1), differenceSwingby2BeforeDsm(:,2),'color', [0.9290, 0.6940, 0.1250]);
plot(differenceSwingby2AfterDsm(:,1), differenceSwingby2AfterDsm(:,2),'color', [0.9290, 0.6940, 0.1250]);
plot(differenceSwingby3BeforeDsm(:,1), differenceSwingby3BeforeDsm(:,2),'color', [0.4940, 0.1840, 0.5560]);
plot(differenceSwingby3AfterDsm(:,1), differenceSwingby3AfterDsm(:,2),'color', [0.4940, 0.1840, 0.5560]);
hold off;
grid;
xlabel('Time [s]');
ylabel('Difference in position [m]');
title('Difference between patched conics and full problem');





% Perturbed case: propagation of the spacecraft within a perturbed environment
clear all;

propagatedStateDepartureBeforeDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_1_leg_0.dat');
propagatedStateDepartureAfterDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_1_leg_1.dat');
propagatedStateSwingby1BeforeDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_1_leg_2.dat');
propagatedStateSwingby1AfterDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_1_leg_3.dat');
propagatedStateSwingby2BeforeDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_1_leg_4.dat');
propagatedStateSwingby2AfterDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_1_leg_5.dat');
propagatedStateSwingby3BeforeDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_1_leg_6.dat');
propagatedStateSwingby3AfterDsm = importdata('../SimulationOutput/fullProblemInterplanetaryTrajectory_1_leg_7.dat');

patchedConicsDepartureBeforeDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_1_leg_0.dat');
patchedConicsDepartureAfterDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_1_leg_1.dat');
patchedConicsSwingby1BeforeDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_1_leg_2.dat');
patchedConicsSwingby1AfterDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_1_leg_3.dat');
patchedConicsSwingby2BeforeDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_1_leg_4.dat');
patchedConicsSwingby2AfterDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_1_leg_5.dat');
patchedConicsSwingby3BeforeDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_1_leg_6.dat');
patchedConicsSwingby3AfterDsm = importdata('../SimulationOutput/patchedConicsInterplanetaryTrajectory_1_leg_7.dat');

differenceDepartureLegBeforeDsm = [propagatedStateDepartureBeforeDsm(:,1) sqrt( (propagatedStateDepartureBeforeDsm(:,2) - patchedConicsDepartureBeforeDsm(:,2)).^2 + (propagatedStateDepartureBeforeDsm(:,3) - patchedConicsDepartureBeforeDsm(:,3)).^2 + (propagatedStateDepartureBeforeDsm(:,4) - patchedConicsDepartureBeforeDsm(:,4)).^2)];
differenceDepartureLegAfterDsm = [propagatedStateDepartureAfterDsm(:,1) sqrt( (propagatedStateDepartureAfterDsm(:,2) - patchedConicsDepartureAfterDsm(:,2)).^2 + (propagatedStateDepartureAfterDsm(:,3) - patchedConicsDepartureAfterDsm(:,3)).^2 + (propagatedStateDepartureAfterDsm(:,4) - patchedConicsDepartureAfterDsm(:,4)).^2)];
differenceSwingby1BeforeDsm = [propagatedStateSwingby1BeforeDsm(:,1) sqrt( (propagatedStateSwingby1BeforeDsm(:,2) - patchedConicsSwingby1BeforeDsm(:,2)).^2 + (propagatedStateSwingby1BeforeDsm(:,3) - patchedConicsSwingby1BeforeDsm(:,3)).^2 + (propagatedStateSwingby1BeforeDsm(:,4) - patchedConicsSwingby1BeforeDsm(:,4)).^2)];
differenceSwingby1AfterDsm = [propagatedStateSwingby1AfterDsm(:,1) sqrt( (propagatedStateSwingby1AfterDsm(:,2) - patchedConicsSwingby1AfterDsm(:,2)).^2 + (propagatedStateSwingby1AfterDsm(:,3) - patchedConicsSwingby1AfterDsm(:,3)).^2 + (propagatedStateSwingby1AfterDsm(:,4) - patchedConicsSwingby1AfterDsm(:,4)).^2)];
differenceSwingby2BeforeDsm = [propagatedStateSwingby2BeforeDsm(:,1) sqrt( (propagatedStateSwingby2BeforeDsm(:,2) - patchedConicsSwingby2BeforeDsm(:,2)).^2 + (propagatedStateSwingby2BeforeDsm(:,3) - patchedConicsSwingby2BeforeDsm(:,3)).^2 + (propagatedStateSwingby2BeforeDsm(:,4) - patchedConicsSwingby2BeforeDsm(:,4)).^2)];
differenceSwingby2AfterDsm = [propagatedStateSwingby2AfterDsm(:,1) sqrt( (propagatedStateSwingby2AfterDsm(:,2) - patchedConicsSwingby2AfterDsm(:,2)).^2 + (propagatedStateSwingby2AfterDsm(:,3) - patchedConicsSwingby2AfterDsm(:,3)).^2 + (propagatedStateSwingby2AfterDsm(:,4) - patchedConicsSwingby2AfterDsm(:,4)).^2)];
differenceSwingby3BeforeDsm = [propagatedStateSwingby3BeforeDsm(:,1) sqrt( (propagatedStateSwingby3BeforeDsm(:,2) - patchedConicsSwingby3BeforeDsm(:,2)).^2 + (propagatedStateSwingby3BeforeDsm(:,3) - patchedConicsSwingby3BeforeDsm(:,3)).^2 + (propagatedStateSwingby3BeforeDsm(:,4) - patchedConicsSwingby3BeforeDsm(:,4)).^2)];
differenceSwingby3AfterDsm = [propagatedStateSwingby3AfterDsm(:,1) sqrt( (propagatedStateSwingby3AfterDsm(:,2) - patchedConicsSwingby3AfterDsm(:,2)).^2 + (propagatedStateSwingby3AfterDsm(:,3) - patchedConicsSwingby3AfterDsm(:,3)).^2 + (propagatedStateSwingby3AfterDsm(:,4) - patchedConicsSwingby3AfterDsm(:,4)).^2)];



figure();
hold on;

% Patched conics trajectory
plot(patchedConicsDepartureBeforeDsm(:,2), patchedConicsDepartureBeforeDsm(:,3), 'color', [0, 0.4470, 0.7410], 'lineWidth', 1.1);
plot(patchedConicsDepartureAfterDsm(:,2), patchedConicsDepartureAfterDsm(:,3), 'color', [0, 0.4470, 0.7410], 'lineWidth', 1.1);
plot(patchedConicsSwingby1BeforeDsm(:,2), patchedConicsSwingby1BeforeDsm(:,3),'color', [0.8500, 0.3250, 0.0980], 'lineWidth', 1.1);
plot(patchedConicsSwingby1AfterDsm(:,2), patchedConicsSwingby1AfterDsm(:,3), 'color', [0.8500, 0.3250, 0.0980], 'lineWidth', 1.1);
plot(patchedConicsSwingby2BeforeDsm(:,2), patchedConicsSwingby2BeforeDsm(:,3), 'color', [0.9290, 0.6940, 0.1250], 'lineWidth', 1.1);
plot(patchedConicsSwingby2AfterDsm(:,2), patchedConicsSwingby2AfterDsm(:,3), 'color', [0.9290, 0.6940, 0.1250], 'lineWidth', 1.1);
plot(patchedConicsSwingby3BeforeDsm(:,2), patchedConicsSwingby3BeforeDsm(:,3), 'color', [0.4940, 0.1840, 0.5560], 'lineWidth', 1.1);
plot(patchedConicsSwingby3AfterDsm(:,2), patchedConicsSwingby3AfterDsm(:,3), 'color', [0.4940, 0.1840, 0.5560], 'lineWidth', 1.1);

% Propagation of the spacecraft 
plot(propagatedStateDepartureBeforeDsm(:,2), propagatedStateDepartureBeforeDsm(:,3), '--', 'color', [0, 0.4470, 0.7410], 'lineWidth', 0.6);
plot(propagatedStateDepartureAfterDsm(:,2), propagatedStateDepartureAfterDsm(:,3), '--', 'color', [0, 0.4470, 0.7410], 'lineWidth', 0.6);
plot(propagatedStateSwingby1BeforeDsm(:,2), propagatedStateSwingby1BeforeDsm(:,3), '--','color', [0.8500, 0.3250, 0.0980], 'lineWidth', 0.6);
plot(propagatedStateSwingby1AfterDsm(:,2), propagatedStateSwingby1AfterDsm(:,3), '--', 'color', [0.8500, 0.3250, 0.0980], 'lineWidth', 0.6);
plot(propagatedStateSwingby2BeforeDsm(:,2), propagatedStateSwingby2BeforeDsm(:,3), '--', 'color', [0.9290, 0.6940, 0.1250], 'lineWidth', 0.6);
plot(propagatedStateSwingby2AfterDsm(:,2), propagatedStateSwingby2AfterDsm(:,3), '--', 'color', [0.9290, 0.6940, 0.1250], 'lineWidth', 0.6);
plot(propagatedStateSwingby3BeforeDsm(:,2), propagatedStateSwingby3BeforeDsm(:,3), '--', 'color', [0.4940, 0.1840, 0.5560], 'lineWidth', 0.6);
plot(propagatedStateSwingby3AfterDsm(:,2), propagatedStateSwingby3AfterDsm(:,3), '--', 'color', [0.4940, 0.1840, 0.5560], 'lineWidth', 0.6);

% Transfer bodies of the trajectory
plot(patchedConicsDepartureBeforeDsm(1,2), patchedConicsDepartureBeforeDsm(1,3), '.k', 'markerSize', 12);
plot(patchedConicsSwingby1BeforeDsm(1,2), patchedConicsSwingby1BeforeDsm(1,3), '.k', 'markerSize', 12);
plot(patchedConicsSwingby2BeforeDsm(1,2), patchedConicsSwingby2BeforeDsm(1,3), '.k', 'markerSize', 12);
plot(patchedConicsSwingby3BeforeDsm(1,2), patchedConicsSwingby3BeforeDsm(1,3), '.k', 'markerSize', 12);
plot(patchedConicsSwingby3AfterDsm(end,2), patchedConicsSwingby3AfterDsm(end,3), '.k', 'markerSize', 12)

% DSM locations along the trajectory
plot(patchedConicsDepartureBeforeDsm(end,2), patchedConicsDepartureBeforeDsm(end,3), '.r', 'markerSize', 12);
plot(patchedConicsSwingby1BeforeDsm(end,2), patchedConicsSwingby1BeforeDsm(end,3), '.r', 'markerSize', 12);
plot(patchedConicsSwingby2BeforeDsm(end,2), patchedConicsSwingby2BeforeDsm(end,3), '.r', 'markerSize', 12);
plot(patchedConicsSwingby3BeforeDsm(end,2), patchedConicsSwingby3BeforeDsm(end,3), '.r', 'markerSize', 12);

hold off;
grid;
xlabel('x [m]', 'fontSize', 12 );
ylabel('y [m]', 'fontSize', 12 );
title('Full problem in the perturbed case');



% Difference between patched conics and propagated state
figure();
semilogy(differenceDepartureLegBeforeDsm(:,1), differenceDepartureLegBeforeDsm(:,2), 'color', [0, 0.4470, 0.7410]);
hold on;
plot(differenceDepartureLegAfterDsm(:,1), differenceDepartureLegAfterDsm(:,2),'color', [0, 0.4470, 0.7410]);
plot(differenceSwingby1BeforeDsm(:,1), differenceSwingby1BeforeDsm(:,2),'color', [0.8500, 0.3250, 0.0980]);
plot(differenceSwingby1AfterDsm(:,1), differenceSwingby1AfterDsm(:,2),'color', [0.8500, 0.3250, 0.0980]);
plot(differenceSwingby2BeforeDsm(:,1), differenceSwingby2BeforeDsm(:,2),'color', [0.9290, 0.6940, 0.1250]);
plot(differenceSwingby2AfterDsm(:,1), differenceSwingby2AfterDsm(:,2),'color', [0.9290, 0.6940, 0.1250]);
plot(differenceSwingby3BeforeDsm(:,1), differenceSwingby3BeforeDsm(:,2),'color', [0.4940, 0.1840, 0.5560]);
plot(differenceSwingby3AfterDsm(:,1), differenceSwingby3AfterDsm(:,2),'color', [0.4940, 0.1840, 0.5560]);
hold off;
grid;
xlabel('Time [s]');
ylabel('Difference in position [m]');
title('Difference between patched conics and full problem');

