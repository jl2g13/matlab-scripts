% function to calculate daylength from lat (in degrees) and year day.

function [daylength] = daylength(lat, yDay)

gamma = double(lat) / double(180) * pi;  % convert lattitude to radians
psi = double(yDay) / double(365) * double(2) * pi; % converting date into an angle


solardec = (0.39637 - 22.9133 .* cos(psi) + 4.02543 .* sin(psi) - 0.38720 .* cos(2*psi) + 0.05200 .* sin(2*psi)) * pi / double(180); % calculating solar declination (Kirk 3rd Ed eq (2.6) p39)

    % calculating daylength
r = -tan(gamma) .* tan(solardec); % see kirk Ed3 eq(2.15)

% alternative to if clause below:
daylength = r * 0;
q1 = find(r <= -1); % q1 in range <= -1
daylength(q1) = 24;
q2 = find(r < 1);   % q2 in range < 1
q3 = setdiff(q2, q1); % selects part of range that is between 1 and -1, i.e. in q2 but not q1
daylength(q3) = (24 / pi) .* acos(r(q3))



% Original function in C: from http://orca.science.oregonstate.edu/faq01.php
% 
% float LatToDayLength(float lat, int yDay)
% {
%  // get lat in radians  % comments in C: use slashes
%  double gamma = lat / 180.0 * M_PI; % double precision for var 'gamma'. 'M_PI' is the irrational number pi in C
% 
%  // convert date into an angle
%  double psi = yDay / 365.0 * 2.0 * M_PI;
% 
%  // calc solar declination
%  // Kirk page 35
%  double solarDec = (0.39637
%             - 22.9133 * cos(psi)
%             + 4.02543 * sin(psi)
%             - 0.38720 * cos(2*psi)
%             + 0.05200 * sin(2*psi)) *M_PI / 180.0;
% 
%  double r = -tan(gamma) * tan(solarDec);
%  if(r <= -1)
%    return 24.0;         % returns the value 24.0 to the calling function ('lattodaylength')
%  else
%    if(fabs(r) < 1) 
%      return 24.0 * acos(r) / M_PI;
%    else
%      return 0;
% }

