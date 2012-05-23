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

% Set simulation data files directory.
simulationDataFilesDirectory = 'build';

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

for i = 1:30
    simulationDataFiles{i} = csvread(strcat(simulationDataFilesDirectory,...
                                            '/galileoSatellite', num2str(i),...
                                            '.dat'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot orbits of Asterix and Obelix.
% Simulation data is given in m and m/s, so needs to be converted to km.

% Plot orbits.
grid on;
xlabel( 'Cartesian x-position [km]' );
ylabel( 'Cartesian y-position [km]' );
zlabel( 'Cartesian z-position [km]' );

for i = 1:30
plot3( simulationDataFiles{i}(:,2)/1000,simulationDataFiles{i}(:,3)/1000,...
    simulationDataFiles{i}(:,4)/1000, 'LineWidth', 3 );
end


for j = 1:size(simulationDataFiles{1},1)
    for i = 1:30
        
        % Plot instantaneous position of Asterix and Obelix.
        figureHandle(i) = plot3( simulationDataFiles{i}(j,2)/1000,...
            simulationDataFiles{i}(j,3)/1000,...
            simulationDataFiles{i}(j,4)/1000, 'o',...
            'MarkerFaceColor', 'y',...
            'MarkerEdgeColor', 'k',...
            'MarkerSize', 10 );
    end
    
    % Capture frame for movie.
    simulationMovieFrames(j) = getframe;
    
    delete( figureHandle );
end

% Play movie.
movie(simulationMovieFrames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%