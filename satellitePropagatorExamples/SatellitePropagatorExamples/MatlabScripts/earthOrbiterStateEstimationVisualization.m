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
observationTimes = (observationTimes - observationTimes(1) )/86400;
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
xlabel('Parameter index [-]')
ylabel('Parameter index [-]')
title('Correlation matrix')

%%%% VISUALIZE WEIGHTED PARTIALS MATRIX: THE HIGHER THE VALUE, THE GREATER THE CONTRUBUTION OF A SINGLE TERM TO THE ESTIMATION
figure(2)
imagesc( log10(weightedPartialsMatrix ));
xlabel('Parameter index [-]')
ylabel('Observation index [-]')
title('log_{10} of (partial derivative matrix scaled by sqrt(weight))')

%%%% PLOT RESIDUALS, SCALED BY SQUARE ROOT OF WEIGHTS, FOR EACH ITERATION OF THE ESTIMATION, COLORED BY LINK ENDS
figure(3)
for i=1:4
    subplot(2,2,i)
    gscatter(observationTimes/86, residualHistory(:,i).*sqrt(weightsDiagonal),observationLinkEnds,'rbk','...xxx',5,'off')
    grid on
    xlabel('Observation time [days]')
    ylabel('Observation residual, scaled by sqrt(weight) [-]' )
    title(strcat('Iteration ',num2str(i)))
    if( i == 2 )
        legend('Range','Angular position','Doppler','Location','NorthWest')
    end
end
suptitle('Observation residual as a function of time, scaled by sqrt(weight), colored by link ends')

%%%% PLOT RESIDUALS, SCALED BY SQUARE ROOT OF WEIGHTS, FOR EACH ITERATION OF THE ESTIMATION, COLORED BY OBSERVABLE TYPE
figure(4)
for i=1:4
    subplot(2,2,i)
    gscatter(observationTimes, residualHistory(:,i).*sqrt(weightsDiagonal),observationTypes,'rbk','...',5,'off')
    grid on
    xlabel('Observation time [days]')
    ylabel('Observation residual, scaled by sqrt(weight) [-]' )
    title(strcat('Iteration ',num2str(i)))
    if( i == 2 )
        legend('Link Ends Set 1','Link Ends Set 2','Link Ends Set 3','Location','NorthWest')
    end
end
suptitle('Observation residual as a function of time, scaled by sqrt(weight), colored by observable type')

%%%% PLOT RESIDUAL HISTOGRAM, SCALED BY SQUARE ROOT OF WEIGHTS, FOR EACH ITERATION OF THE ESTIMATION
figure(5)
for i=1:4
    subplot(2,2,i)
    hist( residualHistory(:,i).*sqrt(weightsDiagonal),40);
    grid on
    ylabel('Number of occurences [-]')
    xlabel('Observation residual, scaled by sqrt(weight) [-]' )
    title(strcat('Iteration ',num2str(i)))
end
suptitle('Observation residual as a function of time, scaled by sqrt(weight)')

%%%% PLOT RATIO OF TRUE TO FORMAL ERROR, AS FUNCTION OF VARIABLE NUMBER,
%%%% AND AS HISTOGRAM. IN IDEAL CASE, DISTRIBUTION IS GAUSSIAN WITH 0 MEAN
%%%% AND 1 STANDARD DEVIATION
figure(6)
subplot(1,2,1)
scatter(1:size(trueError,1),trueError./formalError)
xlabel('Parameter index [-]')
ylabel('True-to-formal error ratio')
grid on

subplot(1,2,2)
hist(trueError./formalError)
ylabel('Number of occurences [-]')
xlabel('True-to-formal error ratio [-]')
grid on

for j=1:6
    set( figure(j), 'Units', 'normalized', 'Position', [0,0,1,1]);
    set( figure(j),'PaperUnits','centimeters','PaperPosition',[0 0 60 40]);
end

sleep(2.0)

for j=1:6
    figure(j)
    if(j==1)
        saveas(gcf,strcat('correlations'),'png');
    elseif(j==2)
        saveas(gcf,strcat('weightedPartials'),'png');
    elseif(j==3)
        saveas(gcf,strcat('residualsPerIterationPerLinkEnd'),'png');
    elseif(j==4)
        saveas(gcf,strcat('residualsPerIterationPerObservable'),'png');
    elseif(j==5)
        saveas(gcf,strcat('residualsPerIterationHistogram'),'png');
    elseif(j==6)
        saveas(gcf,strcat('formalErrorRatios'),'png');
    end
end