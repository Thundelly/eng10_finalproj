fprintf('\n  < Results >\n');

fprintf('\nDeparture planet             ');
disp(planets(depart_num,:));
fprintf('Departure calendar date      ');
disp(calendarDateStr_1);
fprintf('Departure universal time     ');
disp(universalTimeStr_1);
fprintf('\nDeparture julian date        %12.6f', jd1);

fprintf('\n\nArrival planet               ');
disp(planets(arriv_num,:));
fprintf('Arrival calendar date        ');
disp(calendarDateStr_2);
fprintf('Arrival universal time       ');
disp(universalTimeStr_2);
fprintf('\nArrival julian date          %12.6f', jd2');

fprintf('\n\nTransfer time              %12.6f  days \n ', tof);

fprintf('\n\nHeliocentric ecliptic orbital elements of the departure planet\n');
fprintf('--------------------------------------------------------------');
fprintf ('\n        sma (AU)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe1(7)/au, oe1(2), oe1(4), oe1(5));
fprintf ('\n       raan (deg)          true anomaly (deg)       longper(deg)              period (days)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe1(3), oe1(6), oe1(8), T1);

fprintf('\n\nHeliocentric ecliptic orbital elements of the transfer orbit\n')
fprintf('prior to reaching the sphere of influence of the arrival planet\n');
fprintf('---------------------------------------------------------------');
fprintf ('\n        sma (AU)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe(7)/au, oe(2), oe(4)/deg, oe(5)/deg);
fprintf ('\n       raan (deg)          true anomaly (deg)       period (min)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe(3)/deg, oe(6)/deg, T3);

fprintf('\n\nHeliocentric ecliptic orbital elements of the transfer orbit\n')
fprintf('prior to reaching the sphere of influence of the arrival planet\n');
fprintf('---------------------------------------------------------------');
fprintf ('\n        sma (AU)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe3(7)/au, oe3(2), oe3(4)/deg, oe3(5)/deg);
fprintf ('\n       raan (deg)          true anomaly (deg)       period (min)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe3(3)/deg, oe3(6)/deg, T3);

fprintf('\n\nHeliocentric ecliptic orbital elements of the arrival planet\n');
fprintf('--------------------------------------------------------------');
fprintf ('\n        sma (AU)              eccentricity          inclination (deg)         argper (deg)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe2(7)/au, oe2(2), oe2(4), oe2(5));
fprintf ('\n       raan (deg)          true anomaly (deg)       longper(deg)              period (days)');
fprintf ('\n %+16.14e  %+16.14e  %+16.14e  %+16.14e \n', oe2(3), oe2(6), oe2(8), T2);

fprintf('\n\nDeparture velocity vector and magnitude\n');
fprintf('\nx-component of departure velocity                          %12.6f  km/s',V1(1));
fprintf('\ny-component of departure velocity                          %12.6f  km/s', V1(2));
fprintf('\nz-component of departure velocity                          %12.6f  km/s', V1(3));
fprintf('\ndeparture velocity magnitude                               %12.6f  km/s',norm(V1(3)));

fprintf('\n\nArrival velocity vector and magnitude\n');
fprintf('\nx-component of arrival velocity                            %12.6f  km/s', V2(1));
fprintf('\ny-component of arrival velocity                            %12.6f  km/s', V2(2));
fprintf('\nz-component of arrival velocity                            %12.6f  km/s', V2(3));
fprintf('\ndeparture velocity magnitude                               %12.6f  km/s',norm(V2));

fprintf('\n\nHyperbolic excess velocity vector and magnitude at departure\n');
fprintf('\nx-component of hyperbolic excess velocity at departure     %12.6f  km/s',vinf1(1));
fprintf('\ny-component of hyperbolic excess velocity at departure     %12.6f  km/s', vinf1(2));
fprintf('\nz-component of hyperbolic excess velocity at departure     %12.6f  km/s', vinf1(3));
fprintf('\nhyperbolic excess velocity at departure magnitude          %12.6f  km/s',norm(vinf1));

fprintf('\n\nHyperbolic excess velocity vector and magnitude at arrival\n');
fprintf('\nx-component of hyperbolic excess velocity at arrival       %12.6f  km/s', vinf2(1));
fprintf('\ny-component of hyperbolic excess velocity at arrival       %12.6f  km/s', vinf2(2));
fprintf('\nz-component of hyperbolic excess velocity at arrival       %12.6f  km/s', vinf2(3));
fprintf('\nhyperbolic excess velocity at arrival magnitude            %12.6f  km/s',norm(vinf2));

fprintf('\n\n<Planetary departure parameters>\n');

fprintf('\nAltitude of the parking orbit                                             %12.6f  km', a_parking);
fprintf('\nPeriod of the parking orbit                                               %12.6f  min', T_parking);
fprintf('\nSpeed of the space vehicle in its circular orbit                          %12.6f  km/s', vC1);
fprintf('\nRadius to periapsis of the departure hyperbola                            %12.6f  km', rp1);
fprintf('\nEccentricity of the departure hyperbola                                   %12.6f', e_dep);
fprintf('\nSpeed of the space vehicle at the periapsis of the departure hyperbola    %12.6f  km/s', vp1);
fprintf('\nDelta_v for departure                                                     %12.6f  km/s', delta_v_departure);

fprintf('\n\n<Planetary rendezvous parameters>\n');

fprintf('\nAltitude of the capture orbit                                             %12.6f  km', r_capture);
fprintf('\nPeriod of the capture orbit                                               %12.6f  min', T_parking2);
fprintf('\nSpeed of the space vehicle in its circular orbit                          %12.6f  km/s', vC2);
fprintf('\nRadius to periapsis of the arrival hyperbola                              %12.6f  km', r_p_arrival);
fprintf('\nEccentricity of the arrival   hyperbola                                   %12.6f', e_arrive);
fprintf('\nSpeed of the space vehicle at the periapsis of the arrival hyperbola      %12.6f  km/s', vp2);
fprintf('\nDelta_v for arrival                                                       %12.6f  km/s', delta_v_arrival);

fprintf('\n\nTotal delta_v for the mission                                           %12.6f  km/s\n', delta_v_total);