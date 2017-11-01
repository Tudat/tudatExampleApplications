clc
close all
clear all

%% 
%%%%% DEFINE DATA DIRECTORY AND LOAD FILES
dataDirectory = '../SimulationOutput/EarthOrbiterStateEstimationExample/';

correlations = load(strcat(dataDirectory,'earthOrbitEstimationCorrelations.dat'));
residualHistory = load(strcat(dataDirectory,'earthOrbitResidualHistory.dat'));
weightsDiagonal = load(strcat(dataDirectory,'earthOrbitEstimationWeightsDiagonal.dat'));
observations = load(strcat(dataDirectory,'earthOrbitObservationMeasurements.dat'));
observationTimes = load(strcat(dataDirectory,'earthOrbitObservationTimes.dat'));
observationLinkEnds = load(strcat(dataDirectory,'earthOrbitObservationLinkEnds.dat'));
observationTypes = load(strcat(dataDirectory,'earthOrbitObservationObservableTypes.dat'));
partialsMatrix = load(strcat(dataDirectory,'earthOrbitEstimationInformationMatrix.dat'));
trueError = load(strcat(dataDirectory,'earthOrbitObservationTrueEstimationError.dat'));
formalError = load(strcat(dataDirectory,'earthOrbitObservationFormalEstimationError.dat'));

%%%%% CREATE MATRIX OF PARTIAL DERIVATIVES, WITH TERMS SCALED BY SQUARE ROOT OF WEIGHT

weightedPartialsMatrix = partialsMatrix;
for i=1:size(partialsMatrix,2)
    weightedPartialsMatrix(:,i) = ( abs( sqrt( weightsDiagonal ).*weightedPartialsMatrix(:,i ) ) );
end

%%

%%%% VISUALIZE CORRELATION MATRIX
close all
figure(1)
imagesc( abs( correlations ) );

%%%% VISUALIZE WEIGHTED PARTIALS MATRIX: THE HIGHER THE VALUE, THE GREATER THE CONTRUBUTION OF A SINGLE TERM TO THE ESTIMATION
figure(2)
imagesc( weightedPartialsMatrix );

%%%% PLOT RESIDUALS, SCALED BY SQUARE ROOT OF WEIGHTS, FOR EACH ITERATION OF THE ESTIMATION, COLORED BY LINK ENDS
figure(3)
for i=1:4
    subplot(2,2,i)
    gscatter(observationTimes, residualHistory(:,i).*sqrt(weightsDiagonal),observationLinkEnds,'rbk','...xxx',5,'off')
    grid on
end

%%%% PLOT RESIDUALS, SCALED BY SQUARE ROOT OF WEIGHTS, FOR EACH ITERATION OF THE ESTIMATION, COLORED BY OBSERVABLE TYPE
figure(4)
for i=1:4
    subplot(2,2,i)
    gscatter(observationTimes, residualHistory(:,i).*sqrt(weightsDiagonal),observationTypes,'rbk','...',5,'off')
    grid on
end

%%%% PLOT RESIDUAL HISTOGRAM, SCALED BY SQUARE ROOT OF WEIGHTS, FOR EACH ITERATION OF THE ESTIMATION
figure(5)
for i=1:4
     subplot(2,2,i)
    hist( residualHistory(:,i).*sqrt(weightsDiagonal),40);
    grid on
end

%%%% PLOT RATIO OF TRUE TO FORMAL ERROR, AS FUNCTION OF VARIABLE NUMBER,
%%%% AND AS HISTOGRAM. IN IDEAL CASE, DISTRIBUTION IS GAUSSIAN WITH 0 MEAN
%%%% AND 1 STANDARD DEVIATION
figure(6)
subplot(1,2,1)
scatter(1:size(trueError,1),trueError./formalError)

subplot(1,2,2)
hist(trueError./formalError)


