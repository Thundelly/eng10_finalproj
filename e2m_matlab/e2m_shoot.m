function [f, g] = e2m_shoot(x)

% b-plane/orbital element/EI targeting function

% simple shooting method

% input

%  x(1) = current guess for launch v-infinity magnitude (km/sec)
%  x(2) = current guess for launch hyperbola rla (radians)
%  x(3) = current guess for launch hyperbola dla (radians)

% output

%  f(1) = objective function (v-infinity magnitude)
%  f(2) = bdott error
%  f(3) = bdotr error

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global aunit emu pmu rtd rgeo_soi vgeo_soi

global req_mars jdtdb_tip rpt xinct rhyper vhyper

global tevent yevent jdtdb_soi itar_type fpa_tar rmag_tar

global jdtdb_ca rm2sc vm2sc xinc_po rpmag_po

global vinf rla dla rsc_soi vsc_soi thetat

global rsc_ca vsc_ca bdott_user bdotr_user

% current vinf, rla and dla

vinf = x(1);

rla = x(2);

dla = x(3);

% objective function (launch hyperbola v-infinity magnitude)

f(1) = x(1);

% current state vector of launch hyperbola

[rhyper, vhyper] = launch(emu, vinf, dla, rla, rpmag_po, xinc_po);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve for geocentric sphere-of-influence conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up for ode45

options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10, 'Events', @soi_event);

% define maximum search time (seconds)

tof = 25.0 * 86400.0;

[~, ~, tevent, yevent, ~] = ode45(@e2m_eqm1, [0 tof], [rhyper vhyper], options);

jdtdb_soi = jdtdb_tip + tevent / 86400.0;

% spacecraft geocentric state vector at soi (km & km/sec)

rgeo_soi = yevent(1:3);

vgeo_soi = yevent(4:6);

% heliocentric state vector of the spacecraft at soi (au & au/day)

svsun = jplephem(jdtdb_soi, 3, 11);

rsun = svsun(1:3);

vsun = svsun(4:6);

for i = 1:1:3
    
    rsc_soi(i) = rsun(i) + yevent(i) / aunit;
    
    vsc_soi(i) = vsun(i) + 86400.0 * yevent(i + 3) / aunit;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve for closest approach conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up for ode45

options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10, 'Events', @e2m_fpa_event);

% define maximum search time (days)

tof = 500.0;

[~, ~, tevent, yevent, ~] = ode45(@e2m_eqm2, [0 tof], [rsc_soi vsc_soi], options);

rsc_ca = yevent(1:3);

vsc_ca = yevent(4:6);

% julian date at closest approach

jdtdb_ca = jdtdb_soi + tevent;

% state vector of mars at closest approach (au & au/day)

svmars = jplephem(jdtdb_ca, 4, 11);

rmars = svmars(1:3);

vmars = svmars(4:6);

% mars-centered state vector of the spacecraft at closest approach

rm2sc = yevent(1:3)' - rmars(1:3);

vm2sc = yevent(4:6)' - vmars(1:3);

tmatrix = mme2000(jdtdb_ca);

rsc = tmatrix * rm2sc;

vsc = tmatrix * vm2sc;

% mars-to-spacecraft state vector (km and km/sec)

rm2sc = aunit * rsc;

vm2sc = aunit * vsc / 86400.0;

% compute b-plane coordinates

[bplane, ~, ~, ibperr] = rv2bp2(pmu, rm2sc, vm2sc);

% error check

if (ibperr == 1)
    
    snsetStatus(-1);
    
    g = [];
    
    return
    
end

if (itar_type == 1)
    
    % ------------------------------------------
    % target to user-defined B-plane coordinates
    % ------------------------------------------
    
    % b dot t error
    
    f(2) = (bplane(7) - bdott_user) / req_mars;
    
    % b dot r error
    
    f(3) = (bplane(8) - bdotr_user) / req_mars;
    
end

if (itar_type == 2)
    
    % -------------------------------------------------------
    % target to user-defined periapsis radius and inclination
    % -------------------------------------------------------
 
    % extract characteristics of incoming hyperbola
    
    decl_asy = bplane(2);
    
    % rasc_asy = bplane(3);
    
    vinf = bplane(1);
        
    btarg = sqrt(2.0 * pmu * rpt / (vinf * vinf) + rpt * rpt);
    
    % cosine of b-plane angle
    
    ctheta = cos(xinct) / cos(decl_asy);
    
    % check for invalid orbital inclination
    
    tmp = 1.0 - ctheta * ctheta;
    
    if (tmp < 0.0)
        
        fprintf('\nb-plane targeting error!!');
        
        fprintf('\n|inclination| must be > |asymptote declination|');
        
        fprintf('\n\nasymptote declination   %12.6f  degrees', rtd * decl_asy);
        
        quit;
        
    end
    
    stheta = -sqrt(1.0 - ctheta * ctheta);
    
    bdott = btarg * ctheta;
    
    bdotr = btarg * stheta;
    
    % b dot t error
    
    f(2) = (bplane(7) - bdott) / req_mars;
    
    % b dot r error
    
    f(3) = (bplane(8) - bdotr) / req_mars;
    
end

if (itar_type == 3)
 
    % ------------------------------------
    % target to user-defined EI conditions
    % ------------------------------------
    
    ctheta = cos(xinct) / cos(bplane(2));
    
    % check for invalid orbital inclination
    
    tmp = 1.0 - ctheta * ctheta;
    
    if (tmp < 0.0)
        
        fprintf('\nb-plane targeting error!!');
        
        fprintf('\n|inclination| must be > |asymptote declination|');
        
        fprintf('\n\nasymptote declination   %12.6f  degrees', rtd * decl_asy);
        
        quit;
        
    end
    
    stheta = -sqrt(1.0 - ctheta * ctheta);
    
    btarget = cos(fpa_tar) * sqrt(2.0 * pmu * rmag_tar ...
        / (bplane(1) * bplane(1)) + rmag_tar * rmag_tar);
    
    % b dot t error
    
    f(2) = (bplane(7) - btarget * ctheta) / req_mars;
    
    % b dot r error
    
    f(3) = (bplane(8) - btarget * stheta) / req_mars;
    
end

if (itar_type == 4)
    
   % -----------------------
   % target to grazing flyby
   % -----------------------

   vinf = bplane(1);

   btarget = sqrt(2.0 * pmu * req_mars / vinf^2 + req_mars^2);

   bpuser(1) = btarget * sin(thetat);

   bpuser(2) = btarget * cos(thetat);
   
   % b dot t equality constraint
    
   f(2) = (bplane(7) - bpuser(2)) / req_mars;

   % b dot r equality constraint
    
   f(3) = (bplane(8) - bpuser(1)) / req_mars;

end

% transpose objective function/constraint vector

f = f';

% no derivatives

g = [];
