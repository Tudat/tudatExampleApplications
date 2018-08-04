fclose('all'); clear all; close all force; profile off; clc; format long g; rng default;

%% Settings

%...Folder name
outputRepository = '../SimulationOutput/PropagatorTypesComparison';

%...Set file name additions
nameAdditionPropagator = { '_cowell', '_encke', '_kepl', '_equi', '_usm7', '_usm6', '_usmem' };
nameAdditionIntegrator = { '_var', '_const' };

%...Create empty cell arrays for output storage
simulationTimes = cell( length( nameAdditionPropagator ), length( nameAdditionIntegrator ) );
cartesianResults = cell( length( nameAdditionPropagator ), length( nameAdditionIntegrator ) );

%...Set legend labels
legendLabels = { 'Cowell', 'Encke', 'Keplerian', 'Equinoctial', 'USM7', 'USM6', 'USMEM' };

%% Read Output

%...Loop over integrators
for i = 1:length( nameAdditionIntegrator )
    %...Loop over propagators
    for p = 1:length( nameAdditionPropagator )
        %...Set file name
        fileName = fullfile( outputRepository, [ 'cartesian', nameAdditionPropagator{ p }, ...
            nameAdditionIntegrator{ i }, '.dat' ] );
        
        %...Read output
        fileID = fopen( fileName, 'r' );
        temporaryResults = textscan( fileID, repmat( '%f', [ 1, 7 ] ), 'CollectOutput', true, ...
            'Delimiter', ',' );
        
        %...Separate time and Cartesian coordinates
        simulationTimes{ p, i } = temporaryResults{ 1 }( :, 1 );
        cartesianResults{ p, i } = temporaryResults{ 1 }( :, 2:end );
    end
end

%...Clean up
clear i p

%% Read Reference Results

%...Read reference trajectory
fileName = fullfile( outputRepository, 'cartesian_ref_var.dat' );

%...Read output
fileID = fopen( fileName, 'r' );
temporaryResults = textscan( fileID, repmat( '%f', [ 1, 7 ] ), 'CollectOutput', true, 'Delimiter', ',' );

%...Separate time and Cartesian coordinates
referenceTime = temporaryResults{ 1 }( :, 1 );
referenceResults = temporaryResults{ 1 }( :, 2:end );

%...Clean up
clear temporaryResults

%% Read Evaluations And Time

%...Read function evaluations
fileName = fullfile( outputRepository, 'functionEvaluations.dat' );
fileID = fopen( fileName, 'r' );
numberOfFunctionEvaluations = textscan( fileID, '%f', 'CollectOutput', true, 'Delimiter', ',' );
numberOfFunctionEvaluations = numberOfFunctionEvaluations{ 1 };

%...Read computation time
fileName = fullfile( outputRepository, 'propagationTime.dat' );
fileID = fopen( fileName, 'r' );
totalComputationTimes = textscan( fileID, '%f', 'CollectOutput', true, 'Delimiter', ',' );
totalComputationTimes = totalComputationTimes{ 1 };

%...Clean up
clear fileName fileID

%% Interpolate Results

%...Interpolate results to match reference orbit
interpolationTime = referenceTime;
interpolatedCartesianResults = cell( size( cartesianResults ) );

%...Loop over integrators
for i = 1:length( nameAdditionIntegrator )
    %...Loop over propagators
    for p = 1:length( nameAdditionPropagator )
        %...Interpolate
        interpolatedCartesianResults{ p, i } = interp1( simulationTimes{ p, i }, cartesianResults{ p, i }, ...
            interpolationTime, 'spline' );
    end
end

%...Clean up
clear i p

%% Compute RMS Error

%...Define RMS error function
rmsErrorFunction = @( result ) [ sqrt( sum( rms( referenceResults( :, 1:3 ) - result( :, 1:3 ) ).^2, 2 ) ), ...
    sqrt( sum( rms( referenceResults( :, 4:6 ) - result( :, 4:6 ) ).^2, 2 ) ) ];

%...Compute RMS errors in position and velocity
rmsErrors = cellfun( rmsErrorFunction, interpolatedCartesianResults, 'UniformOutput', false );
positionRMSErrorVariableStep = cellfun( @(x) x(1), rmsErrors( :, 1 ) );
positionRMSErrorConstantStep = cellfun( @(x) x(1), rmsErrors( :, 2 ) );

%...Plot RMS errors
F = figure( 'rend', 'painters', 'pos', [ 440, 378, 840, 630 ] );
styles = { 'o', 'd', 's', 'v', 'h', 'p', '^' };

subplot( 1, 2, 1 )
hold on
for p = 1:length( nameAdditionPropagator )
    scatter( numberOfFunctionEvaluations( p ), positionRMSErrorVariableStep( p ), 250, 'filled', styles{ p } )
end
hold off
xlabel( 'Number of Function Evaluations [-]' )
ylabel( 'RMS Error in Position [m]' )
grid on
title( 'Variable Step Size' )
set( gca, 'FontSize', 15, 'YScale', 'log' )
% legend( legendLabels{ : }, 'Position', 'Best' )
[ ~, icons ] = legend( legendLabels{:}, 'Location', 'Best' );
for p = 1:length( nameAdditionPropagator )
    icons( p ).FontSize = 12.5;
    icons( p + length( nameAdditionPropagator ) ).Children.MarkerSize = 10;
end

subplot( 1, 2, 2 )
hold on
for p = 1:length( nameAdditionPropagator )
    scatter( totalComputationTimes( p ), positionRMSErrorConstantStep( p ), 250, 'filled', styles{ p } )
end
hold off
xlabel( 'Computation Time [s]' )
ylabel( 'RMS Error in Position [m]' )
grid on
title( 'Constant Step Size' )
set( gca, 'FontSize', 15, 'YScale', 'log' )
% legend( legendLabels{ : }, 'Position', 'Best' )
[ ~, icons ] = legend( legendLabels{:}, 'Location', 'Best' );
for p = 1:length( nameAdditionPropagator )
    icons( p ).FontSize = 12.5;
    icons( p + length( nameAdditionPropagator ) ).Children.MarkerSize = 10;
end

saveas( F, 'propagatorTypesComparison', 'epsc' )

%...Clean up
clear F styles icons p