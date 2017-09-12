clc
close all
clear all

dataDirectory = '../SimulationOutput/EarthOrbiterStateEstimationExample/';

correlations = load(strcat(dataDirectory,'earthOrbitBasicEstimationCorrelations.dat'));
residualHistory = load(strcat(dataDirectory,'earthOrbitBasicResidualHistory.dat'));
weightsDiagonal = load(strcat(dataDirectory,'earthOrbitBasicEstimationWeightsDiagonal.dat'));
observations = load(strcat(dataDirectory,'earthOrbitBasicObservationMeasurements.dat'));
observationTimes = load(strcat(dataDirectory,'earthOrbitBasicObservationTimes.dat'));
observationLinkEnds = load(strcat(dataDirectory,'earthOrbitBasicObservationLinkEnds.dat'));
observationTypes = load(strcat(dataDirectory,'earthOrbitBasicObservationObservableTypes.dat'));
partialsMatrix = load(strcat(dataDirectory,'earthOrbitBasicEstimationInformationMatrix.dat'));
trueError = load(strcat(dataDirectory,'earthOrbitBasicObservationTrueEstimationError.dat'));
formalError = load(strcat(dataDirectory,'earthOrbitBasicObservationFormalEstimationError.dat'));


weightedPartialsMatrix = partialsMatrix;

for i=1:size(partialsMatrix,2)
    weightedPartialsMatrix(:,i) = ( abs( sqrt( weightsDiagonal ).*weightedPartialsMatrix(:,i ) ) );
end

%%
close all
figure(1)
imagesc( abs( correlations ) );

figure(2)
imagesc( weightedPartialsMatrix );

%%
figure(3)
for i=1:4
    subplot(2,2,i)
    gscatter(observationTimes, residualHistory(:,i).*sqrt(weightsDiagonal),observationLinkEnds,'rbk','...xxx',5,'off')
    grid on
end

figure(4)
for i=1:4
    subplot(2,2,i)
    gscatter(observationTimes, residualHistory(:,i).*sqrt(weightsDiagonal),observationTypes,'rbk','...',5,'off')
    grid on
end

figure(5)
for i=1:4
        subplot(2,2,i)
    hist( residualHistory(:,i).*sqrt(weightsDiagonal),40);
    grid on
end

figure(6)
subplot(1,2,1)
scatter(1:size(trueError,1),trueError./formalError)

subplot(1,2,2)
hist(trueError./formalError)


