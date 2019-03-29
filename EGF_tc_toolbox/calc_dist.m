function dist = calc_dist(lat1d,lon1d,lat2d,lon2d)
% Calculates distance in km based on the Haversine formulas
% (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)
% 
% Input: 
%       lat1d = latitude in degrees for station 1
%       lon1d = longitude in degrees for station 1
%       lat2d = latitude in degrees for station 2
%       lon2d = longitude in degrees for station 2
% 
% Output:
%       dist = distance in km between the two stations
%

radius = 6371; % km 

% Find the coodinates in radians 
lat1 = lat1d*pi/180;
lat2 = lat2d*pi/180;
lon1 = lon1d*pi/180;
lon2 = lon2d*pi/180;

deltaLat = lat2-lat1;
deltaLon = lon2-lon1;

a = sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2;
c = 2*atan2(sqrt(a),sqrt(1-a));

dist = radius*c;    % Haversine distance

end