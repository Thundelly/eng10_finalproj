function ydot = e2m_eqm1 (t, y)

% geocentric equations of motion

% includes j2 earth gravity, and point-mass
% gravity of the sun and moon

% input

%  t = simulation time since hyperbolic injection (seconds)
%  y = spacecraft state vector (km & km/sec)

% output

%  ydot = integration vector

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global emu mmu smu req j2 aunit jdtdb_tip

% J2 acceleration due to the earth

r2 = y(1) * y(1) + y(2) * y(2) + y(3) * y(3);

r1 = sqrt(r2);

r3 = r2 * r1;
   
r5 = r2 * r3;
   
d1 = -1.5 * j2 * req * req * emu / r5;
   
d2 = 1.0 - 5.0 * y(3) * y(3) / r2;

agrav(1) = y(1) * (d1 * d2 - emu / r3);
   
agrav(2) = y(2) * (d1 * d2 - emu / r3);
   
agrav(3) = y(3) * (d1 * (d2 + 2) - emu / r3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% point-mass acceleration due to the moon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jdate = jdtdb_tip + t / 86400.0d0;

svmoon = jplephem(jdate, 10, 3);

rmoon = aunit * svmoon(1:3);
   
% compute selenocentric state vector of the spacecraft

rm2sc = y(1:3) - rmoon;

% f(q) formulation

qmoon = dot(y(1:3), y(1:3) - 2.0 * rmoon) / dot(rmoon, rmoon);

fmoon = qmoon * ((3.0 + 3.0 * qmoon + qmoon * qmoon) / (1.0 + (1.0 + qmoon)^1.5));

d3 = norm(rm2sc) * dot(rm2sc, rm2sc);

amoon = -mmu * (y(1:3) + fmoon * rmoon) / d3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% point-mass acceleration due to the sun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

svsun = jplephem(jdate, 11, 3);

rsun = aunit * svsun(1:3);

% compute heliocentric state vector of the spacecraft

rs2sc = y(1:3) - rsun;

% f(q) formulation

qsun = dot(y(1:3), y(1:3) - 2.0 * rsun) / dot(rsun, rsun);

fsun = qsun * ((3.0 + 3.0 * qsun + qsun * qsun) / (1.0 + (1.0 + qsun)^1.5));

d3 = norm(rs2sc) * dot(rs2sc, rs2sc);

asun = -smu * (y(1:3) + fsun * rsun) / d3;

% compute total integration vector

ydot = [ y(4)
         y(5)
         y(6)
         agrav(1) + amoon(1) + asun(1)
         agrav(2) + amoon(2) + asun(2)
         agrav(3) + amoon(3) + asun(3)];
