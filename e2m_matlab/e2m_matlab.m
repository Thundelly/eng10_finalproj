% e2m_matlab_64bit.m     April 16, 2018

% Earth-to-Mars trajectory optimization

% JPL DE421 ephemeris 

% 64 bit SNOPT algorithm (June 17, 2015 version)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

oev2 = zeros(7, 1);

deltav = zeros(3, 1);

global aunit iephem km ephname dtr rtd

global emu smu mmu pmu xmu itar_type rkcoef

global j2 req req_mars rsoi 

global jdtdb0 jdtdb_tip jdtdb_soi jdtdb_ca

global otype rito vito rm2sc vm2sc rfto vfto

global ip1 ip2 dv1 dv2 rpt xinct fpa_tar rmag_tar

global xinc_po rpmag_po rhyper vhyper thetat

global vinf rla dla rsc_soi vsc_soi rsc_ca vsc_ca

global bdott_user bdotr_user rgeo_soi vgeo_soi

% read leap seconds data file

readleap;

% angular conversion factors

rtd = 180.0 / pi;

dtr = pi / 180.0;

atr = dtr / 3600.0;

% de421 value for astronomical unit (kilometers)

aunit = 149597870.699626200;

% gravitational constant of the Earth (km^3/sec^2)

emu = 398600.4415;

% equatorial radius of the Earth (kilometers)

req = 6378.1378;

% j2 gravity coefficient of the Earth

j2 = 0.00108263;

% sphere-of-influence of the Earth (kilometers)

rsoi = 925000.0;

% gravitational constant of the moon (km^3/sec^2)

mmu = 4902.800238;

% gravitational constant of mars (km^3/sec^2)

pmu = 42828.376212;
      
% radius of mars (kilometers)

req_mars = 3396.2;

% initialize de421 ephemeris

ephname = 'de421.bin';

iephem = 1;

km = 0;

% initialize RKF7(8) integration method

rkcoef = 1;

% define "reference" julian date (1/1/2000)

jdtdb0 = 2451544.5;

%*********************************************
% DE421 gravitational constants (au**3/day**2)
%*********************************************

% earth-moon mass ratio

emrat = 0.813005674615454410e+02;

% sun

xmu(1) = 0.295912208285591149e-03;

smu = xmu(1) * aunit^3 / 86400.0^2;

% mercury

xmu(2) = 0.491254745145081187e-10;

% venus

xmu(3) = 0.724345233302178764e-09;

% earth

xmu(4) = 0.899701152970881167e-09;

xmu(4) = emu * 86400.0^2 / aunit^3;

% mars

xmu(5) = 0.954954891857520080e-10;

% jupiter

xmu(6) = 0.282534585557186464e-06;

% saturn

xmu(7) = 0.845972727826697794e-07;

% uranus

xmu(8) = 0.129202491678196939e-07;

% neptune

xmu(9) = 0.152435890078427628e-07;

% read simulation definition data file

[filename, pathname] = uigetfile('*.in', 'Please select the input file to read');

[fid, otype, jdtdb_tip, ddays1, jdtdb_arrival, ddays2, alt_po, aziml, ...
    xlatgs, rpt, xinct, fpa_tar, alt_tar, itar_type, ...
    bdott_user, bdotr_user, thetat] = e2m_readdata(filename);

% set target radius at Mars (kilometers)

rmag_tar = req_mars + alt_tar;

% begin simulation

clc; home;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% step 1 - solve two-body lambert problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% departure and arrival planets

ip1 = 3;

ip2 = 4;

% posigrade transfer

direct = 1;

% less than one rev transfer

revmax = 0;

mu = smu;

if (otype < 4)
    
   %%%%%%%%%%%%%%%%%%%%%%%%%
   % find optimal solution %
   %%%%%%%%%%%%%%%%%%%%%%%%%
      
   xg(1) = jdtdb_tip - jdtdb0;
       
   xg(2) = jdtdb_arrival - jdtdb0;
       
   xg = xg';

   % bounds on control variables

   xlwr(1) = xg(1) - ddays1;
   
   xupr(1) = xg(1) + ddays1;

   xlwr(2) = xg(2) - ddays2;
   
   xupr(2) = xg(2) + ddays2;

   xlwr = xlwr';
   
   xupr = xupr';

   % bounds on objective function

   flow(1) = 0.0;
   
   fupp(1) = +Inf;

   xmul1 = zeros(2, 1);

   xstate1 = zeros(2, 1);

   fmul1 = zeros(1, 1);

   fstate1 = zeros(1, 1);

   % find optimum
       
   snscreen on;

   [x, ~, ~, ~, ~] = snopt(xg, xlwr, xupr, xmul1, xstate1, ...
        flow, fupp, fmul1, fstate1, 'e2m_deltav');

   jdtdb_tip = x(1) + jdtdb0;

   jdtdb_arrival = x(2) + jdtdb0;

   % transfer time (days)
   
   taud = jdtdb_arrival - jdtdb_tip;
  
   % evaluate current solution
   
   [~, ~] = e2m_deltav(x);
   
   dvm1 = norm(dv1);
   
else
    
   %%%%%%%%%%%%%%%%%
   % no optimization
   %%%%%%%%%%%%%%%%%
   
   x(1) = jdtdb_tip - jdtdb0;
       
   x(2) = jdtdb_arrival - jdtdb0;
     
   [~, ~] = e2m_deltav(x);
     
   taud = jdtdb_arrival - jdtdb_tip;
   
   dvm1 = norm(dv1);
   
end

%%%%%%%%%%%%%%%%%
% print results %
%%%%%%%%%%%%%%%%%

% convert solution julian dates to calendar dates and utc

jdutc1 = tdb2utc (jdtdb_tip);

[cdstr1, utstr1] = jd2str(jdutc1);

fprintf('\nEarth-to-Mars mission design\n');

fprintf('\n=========================');
fprintf('\ntwo-body Lambert solution');
fprintf('\n=========================');

switch otype
    
case 1  
    
  fprintf('\n\nminimize departure delta-v\n');
  
case 2
    
  fprintf('\n\nminimize arrival delta-v\n');
  
case 3 
    
  fprintf('\n\nminimize total delta-v\n');
  
case 4  
    
  fprintf('\n\nno optimization\n');
  
end

fprintf('\ndeparture heliocentric delta-v vector and magnitude');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * dv1(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * dv1(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * dv1(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second\n\n', 1000.0 * norm(dv1));

fprintf('\narrival heliocentric delta-v vector and magnitude');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * dv2(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * dv2(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * dv2(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second\n', 1000.0 * norm(dv2));

% print orbital elements of the earth at departure

fprintf('\n\nheliocentric coordinates of the Earth at departure');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nUTC calendar date    ');

disp(cdstr1);

fprintf('\nUTC time             ');

disp(utstr1);

fprintf('\nUTC Julian Date      %12.8f\n', jdutc1);

svearth = jplephem(jdtdb_tip, ip1, 11);

r = svearth(1:3);

v = svearth(4:6);

oev = eci2orb1(smu, aunit * r, aunit * v / 86400.0);

oeprint2(smu, oev);

svprint(aunit * r, aunit * v / 86400.0);

% print orbital elements of the heliocentric transfer orbit

fprintf('\nspacecraft heliocentric coordinates after the first impulse');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nUTC calendar date    ');

disp(cdstr1);

fprintf('\nUTC time             ');

disp(utstr1);

fprintf('\nUTC Julian Date      %12.8f\n', jdutc1);

oev = eci2orb1(smu, rito, vito');

oeprint2(smu, oev);
   
svprint(rito, vito');

fprintf('\nspacecraft heliocentric coordinates prior to the second impulse');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

jdutc2 = tdb2utc (jdtdb_arrival);

[cdstr2, utstr2] = jd2str(jdutc2);

fprintf('\nUTC calendar date    ');

disp(cdstr2);

fprintf('\nUTC time             ');

disp(utstr2);

fprintf('\nUTC Julian Date      %12.8f\n', jdutc2');

% propagate transfer orbit

[rfto, vfto] = twobody2 (smu, 86400.0 * taud, rito, vito);

oev = eci2orb1(smu, rfto, vfto);

oeprint2(smu, oev);
   
svprint(rfto, vfto);

fprintf('\nspacecraft heliocentric coordinates after the second impulse');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nUTC calendar date    ');

disp(cdstr2);

fprintf('\nUTC time             ');

disp(utstr2);

fprintf('\nUTC Julian Date      %12.8f\n', jdutc2');

rscf = rfto;

vscf = vfto + dv2;

oev = eci2orb1(smu, rfto, vfto);

oeprint2(smu, oev);
   
svprint(rscf, vscf);

% print orbital elements of Mars at arrival

fprintf('\nheliocentric coordinates of Mars at arrival');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nUTC calendar date    ');

disp(cdstr2);

fprintf('\nUTC time             ');

disp(utstr2);

fprintf('\nUTC Julian Date      %12.8f\n', jdutc2');

svmars = jplephem(jdtdb_arrival, ip2, 11);

r = svmars(1:3);

v = svmars(4:6);

oev = eci2orb1(smu, aunit * r, aunit * v / 86400.0);

oeprint2(smu, oev);

svprint(aunit * r, aunit * v / 86400.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute orientation of the departure hyperbola
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

decl1 = 90.0 - rtd * acos(dv1(3) / dvm1);

rasc1 = rtd * atan3(dv1(2), dv1(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2 - compute characteristics of launch hyperbola
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rpmag_po = req + alt_po;

xinc_po = acos(cos(xlatgs) * sin(aziml));

% check for invalid orbital inclination
     
if (abs(xinc_po) < abs(dtr * decl1))
    
   fprintf('\npark orbit error!!\n');
       
   fprintf('\n|inclination| must be > |asymptote declination|');
   
   fprintf('\n\npark orbit inclination   %12.6f  degrees', rtd * xinc_po);
   
   fprintf('\n\nasymptote declination    %12.6f  degrees\n\n', decl1);
       
   quit;
   
end
    
[rhyper, vhyper] = launch(emu, dvm1, dtr * decl1, dtr * rasc1, rpmag_po, xinc_po);

oev1 = eci2orb1 (emu, rhyper, vhyper);

% compute orbital elements of circular park orbit

for i = 1:1:6
    
    oev2(i) = oev1(i);
    
end

oev2(1) = rpmag_po;

oev2(2) = 0.0;

oev2(4) = 0.0;

oev2(6) = oev1(4);

[rpo, vpo] = orb2eci (emu, oev2);

% compute injection delta-v

for i = 1:1:3
    
    deltav(i) = vhyper(i) - vpo(i);
    
end

decl1 = 90.0 - rtd * acos(dv1(3) / norm(dv1));

rasc1 = rtd * atan3(dv1(2), dv1(1));

fprintf('\npark orbit and departure hyperbola characteristics');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\npark orbit');
fprintf('\n----------\n');

oeprint1(emu, oev2);

svprint (rpo, vpo);

fprintf('\ndeparture hyperbola');
fprintf('\n-------------------\n');

fprintf('\nc3                          %12.6f  kilometer^2/second^2\n', dvm1 * dvm1);

fprintf('\nv-infinity                  %12.6f  meters/second\n', 1000.0 * dvm1);

fprintf('\nasymptote right ascension   %12.6f  degrees\n', rasc1);

fprintf('\nasymptote declination       %12.6f  degrees\n', decl1);

fprintf('\nperigee altitude            %12.6f  kilometers\n', alt_po);

fprintf('\nlaunch azimuth              %12.6f  degrees\n', rtd * aziml);

fprintf('\nlaunch site latitude        %12.6f  degrees\n\n', rtd * xlatgs);

fprintf('\nUTC calendar date    ');

disp(cdstr1);

fprintf('\nUTC time             ');

disp(utstr1);

fprintf('\nUTC Julian Date      %12.8f\n', jdutc1);

oeprint1(emu, oev1);

svprint (rhyper, vhyper);

fprintf('\nhyperbolic injection delta-v vector and magnitude');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * deltav(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * deltav(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * deltav(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second', 1000.0 * norm(deltav));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve b-plane targeting problem using a simple shooting method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\nplease wait, solving b-plane targeting problem ...\n\n');

% initial guess for launch vinf, rla and dla

xg(1) = dvm1;

xg(2) = dtr * rasc1;

xg(3) = dtr * decl1;

if (otype == 4)
    
   % transpose
   
   xg = xg';
   
end

% define lower and upper bounds for vinf, rla and dla

xlwr(1) = xg(1) - 0.05;
xupr(1) = xg(1) + 0.05;

xlwr(2) = xg(2) - 10.0 * dtr;
xupr(2) = xg(2) + 10.0 * dtr;

xlwr(3) = xg(3) - 1.0 * dtr;
xupr(3) = xg(3) + 1.0 * dtr;

if (otype == 4)
    
   % transpose bounds
   
   xlwr = xlwr';
   
   xupr = xupr';
   
end

% bounds on objective function

flow(1) = 0.0;

fupp(1) = +Inf;

% bounds on final b-plane/orbital element equality constraints

flow(2) = 0.0;
fupp(2) = 0.0;

flow(3) = 0.0;
fupp(3) = 0.0;

flow = flow';

fupp = fupp';

snscreen on;

xmul = zeros(3, 1);

xstate = zeros(3, 1);

fmul = zeros(3, 1);

fstate = zeros(3, 1);
   
[x, ~, inform, xmul, fmul] = snopt(xg, xlwr, xupr, xmul, xstate, ...
        flow, fupp, fmul, fstate, 'e2m_shoot');

% evaluate final solution

[f, g] = e2m_shoot(x);

% compute orbital elements of launch hyperbola

oev1 = eci2orb1 (emu, rhyper, vhyper);

% compute orbital elements of circular park orbit

for i = 1:1:6
    
    oev2(i) = oev1(i);
    
end

oev2(1) = rpmag_po;

oev2(2) = 0.0;

oev2(4) = 0.0;

oev2(6) = oev1(4);

[rpo, vpo] = orb2eci (emu, oev2);

% compute injection delta-v (km/sec)

for i = 1:1:3
    
    deltav(i) = vhyper(i) - vpo(i);
    
end

decl1 = 90.0 - rtd * acos(deltav(3) / norm(deltav));

rasc1 = rtd * atan3(deltav(2), deltav(1));

fprintf('\n=======================');
fprintf('\noptimal n-body solution');
fprintf('\n=======================\n');

if (itar_type == 1)
    
   fprintf('\nB-plane targeting\n');
   
end

if (itar_type == 2)
    
   fprintf('\norbital element targeting\n');
   
end

if (itar_type == 3)
    
   fprintf('\nEI conditions targeting\n');
   
end

if (itar_type == 4)
    
   fprintf('\ngrazing flyby targeting\n');
   
end

fprintf('\npark orbit and departure hyperbola characteristics');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\npark orbit');
fprintf('\n----------\n');

oeprint1(emu, oev2);

svprint (rpo, vpo);

fprintf('\ndeparture hyperbola');
fprintf('\n-------------------\n');

vinf = f(1);

fprintf('\nc3                          %12.6f  km^2/sec^2\n', vinf * vinf);

fprintf('\nv-infinity                  %12.6f  meters/second\n', 1000.0 * vinf);

fprintf('\nasymptote right ascension   %12.6f  degrees\n', rtd * rla);

fprintf('\nasymptote declination       %12.6f  degrees\n\n', rtd * dla);

fprintf('\nUTC calendar date    ');

disp(cdstr1);

fprintf('\nUTC time             ');

disp(utstr1);

fprintf('\nUTC Julian Date      %12.8f\n', jdtdb_tip);

oeprint1(emu, oev1);

svprint (rhyper, vhyper);

fprintf('\nhyperbolic injection delta-v vector and magnitude');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * deltav(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * deltav(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * deltav(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second\n', 1000.0 * norm(deltav));

fprintf('\ntransfer time               %12.6f  days\n', jdtdb_ca - jdtdb_tip);

fprintf('\ntime and conditions at Mars closest approach');
fprintf('\n(Mars mean equator and IAU node of epoch)');
fprintf('\n-----------------------------------------\n');

jdutc_ca = tdb2utc(jdtdb_ca);

[cdstr_ca, utstr_ca] = jd2str(jdutc_ca);

fprintf('\nUTC calendar date    ');

disp(cdstr_ca);

fprintf('\nUTC time             ');

disp(utstr_ca);

fprintf('\nUTC Julian Date      %12.8f\n', jdutc_ca);

% orbital elements at closest approach

oev = eci2orb1(pmu, rm2sc, vm2sc);

oeprint1(pmu, oev);

svprint(rm2sc, vm2sc);

% compute b-plane coordinates

[bplane, rv, tv, ibperr] = rv2bp2(pmu, rm2sc, vm2sc);

% print results

fprintf('\nB-plane coordinates at Mars closest approach');
fprintf('\n(Mars mean equator and IAU node of epoch)');
fprintf('\n-----------------------------------------\n');

fprintf ('\nb-magnitude                %12.6f  kilometers', sqrt(bplane(7)^2 + bplane(8)^2));

fprintf ('\nb dot r                    %12.6f  kilometers', bplane(8));

fprintf ('\nb dot t                    %12.6f  kilometers', bplane(7));

fprintf ('\nb-plane angle              %12.6f  degrees', rtd * bplane(6));

fprintf ('\nv-infinity                 %12.6f  meters/second', 1000.0 * bplane(1));
      
fprintf ('\nr-periapsis                %12.6f  kilometers', bplane(5));
      
fprintf ('\ndecl-asymptote             %12.6f  degrees', rtd * bplane(2));
                  
fprintf ('\nrasc-asymptote             %12.6f  degrees\n', rtd * bplane(3));

% sine of flight path angle at closest approach

sfpa = rm2sc' * vm2sc /(norm(rm2sc) * norm(vm2sc));

fprintf ('\nflight path angle          %12.6f  degrees\n', rtd * asin(sfpa));

fprintf('\nspacecraft heliocentric coordinates at closest approach');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nUTC calendar date    ');

disp(cdstr_ca);

fprintf('\nUTC time             ');

disp(utstr_ca);

fprintf('\nUTC Julian Date      %12.8f\n', jdutc_ca);

oev_ca = eci2orb1 (smu, aunit * rsc_ca, aunit * vsc_ca / 86400.0);

oeprint1(smu, oev_ca);

svprint(aunit * rsc_ca, aunit * vsc_ca / 86400.0);

fprintf('\nheliocentric coordinates of Mars at closest approach');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nUTC calendar date    ');

disp(cdstr_ca);

fprintf('\nUTC time             ');

disp(utstr_ca);

fprintf('\nUTC Julian Date      %12.8f\n', jdutc_ca);

svmars = jplephem(jdtdb_ca, ip2, 11);

r = svmars(1:3);

v = svmars(4:6);

oev = eci2orb1(smu, aunit * r, aunit * v / 86400.0);

oeprint2(smu, oev);

svprint(aunit * r, aunit * v / 86400.0);

fprintf('\nspacecraft geocentric coordinates at the Earth SOI');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

jdutc_soi = tdb2utc(jdtdb_soi);

[cdstr_soi, utstr_soi] = jd2str(jdutc_soi);

fprintf('\nUTC calendar date    ');

disp(cdstr_soi);

fprintf('\nUTC time             ');

disp(utstr_soi);

fprintf('\nUTC Julian Date      %12.8f\n', jdutc_soi);

oev_soi = eci2orb1 (emu, rgeo_soi, vgeo_soi);

oeprint1(emu, oev_soi);

svprint(rgeo_soi, vgeo_soi);

fprintf('\nspacecraft heliocentric coordinates at the Earth SOI');
fprintf('\n(Earth mean equator and equinox of J2000)');
fprintf('\n-----------------------------------------\n');

fprintf('\nUTC calendar date    ');

disp(cdstr_soi);

fprintf('\nUTC time             ');

disp(utstr_soi);

fprintf('\nUTC Julian Date      %12.8f\n', jdutc_soi);

oev_soi = eci2orb1 (smu, aunit * rsc_soi, aunit * vsc_soi / 86400.0);

oeprint1(smu, oev_soi);

svprint(aunit * rsc_soi, aunit * vsc_soi / 86400.0);

while(1)
    
   fprintf('\n\nwould you like to create trajectory graphics (y = yes, n = no)\n');

   slct = input('? ', 's');
   
   if (slct == 'y' || slct == 'n')
       
      break;
      
   end 
   
end

if (slct == 'y')  
    
   % create heliocentric and areocentric trajectory graphics

   marsplot;
   
end
