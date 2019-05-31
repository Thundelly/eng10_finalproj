function [value, isterminal, direction] = soi_event(t, y)

% sphere-of-influence event function

% input

%  t = current simulation time
%  y = current spacecraft geocentric state vector (km & km/sec)

% output

%  value = difference between current position and soi (kilometers)

% global

%  rsoi = Earth sphere-of-influence (kilometers)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global rsoi

% difference between current geocentric distance and soi (kilometers)

value = norm(y(1:3)) - rsoi;

isterminal = 1;

direction =  [];

