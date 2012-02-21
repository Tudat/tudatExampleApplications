%   Copyright (c) 2010-2012 Delft University of Technology.
%
%    This software is protected by national and international copyright.
%    Any unauthorized use, reproduction or modification is unlawful and
%    will be prosecuted. Commercial and non-private application of the
%    software in any form is strictly prohibited unless otherwise granted
%    by the authors.
%
%    The code is provided without any warranty; without even the implied
%    warranty of merchantibility or fitness for a particular purpose.
%
%    Changelog
%      YYMMDD    Author            Comment
%      101111    K. Kumar          File created.
%
%    References
%

clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input deck.

% Set simulation data files.

% Set Asterix simulation data file
asterixSimulationDataFile = 
         
% Set Asterix simulation data file
obelixSimulationDataFile =

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Earth.
% This code was provided by J. Melman, TU Delft.
% Axes are given in km.

% Load data file.
load( 'topo.mat' )

% Set radius of Earth in km.
r_E = 6378.137;

% Create geometric sphere.
[x y z] = sphere( 50 );

% Set plot properties.
hold on
colormap( topomap1 );

props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor = 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;

% Plot Earth.
surface( r_E * x, r_E * y, r_E * z, props );
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read and store simulation data files.
% First column is epoch, subsequent columns are Cartesian state elements.

% Read and store Asterix simulation data
asterixSimulationData = csvread( asterixSimulationDataFile );

% Read and store Obelix simulation data.
obelixSimulationData = csvread( obelixSimulationDataFile );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot orbits of Asterix and Obelix.
% Simulation data is given in m and m/s, so needs to be converted to km.

% Convert data to km.
asterixSimulationData = asterixSimulationData / 1000;
obelixSimulationData = obelixSimulationData / 1000;

% Plot orbits.
grid on;
xlabel( 'Cartesian x-position [km]' );
ylabel( 'Cartesian y-position [km]' );
zlabel( 'Cartesian z-position [km]' );

plot3( asterixSimulationData(:,2),asterixSimulationData(:,3),...
    asterixSimulationData(:,4), 'LineWidth', 3 );
plot3( obelixSimulationData(:,2),obelixSimulationData(:,3),...
    obelixSimulationData(:,4), 'r', 'LineWidth', 3 );

for i = 1:size(asterixSimulationData,1)
    
    % Plot instantaneous position of Asterix and Obelix.
    figureHandleAsterix = plot3( asterixSimulationData(i,2),...
                                 asterixSimulationData(i,3),...
                                 asterixSimulationData(i,4), 'o',...
                                 'MarkerFaceColor', 'y',...
                                 'MarkerEdgeColor', 'k',...
                                 'MarkerSize', 10 );
                             
    figureHandleObelix = plot3( obelixSimulationData(i,2),...
                                obelixSimulationData(i,3),...
                                obelixSimulationData(i,4), 'o',...
                                'MarkerFaceColor', 'y',...
                                'MarkerEdgeColor', 'k',...
                                'MarkerSize', 10 );
    
   % Capture frame for movie.
    simulationMovieFrames(i) = getframe;
    
    delete( figureHandleAsterix );
    delete( figureHandleObelix );
end

% Play movie.
movie(simulationMovieFrames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%