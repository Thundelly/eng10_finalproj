function [rhyper, vhyper] = launch(mu, vinf, decl_asy, rasc_asy, rpmag, xinc)

% orbital elements of a launch hyperbola

% input

%  mu       = gravitational constant (km**3/sec**2)
%  vinf     = v-infinity magnitude (kilometers/second)
%  decl_asy = declination of outgoing asymptote (radians)
%  rasc_asy = right ascension of outgoing asymptote (radians)
%  rpmag    = perigee radius of launch hyperbola (kilometers)
%  xinc     = launch hyperbola inclination (radians)

% output

%  rhyper = position vector at periapsis of launch hyperbola (km)
%  vhyper = velocity vector at periapsis of launch hyperbola (km/sec)

% NOTE: |xinc| must be > |decl_asy|

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pv = [0.0 0.0 1.0];

theta = acos(cos(xinc) / cos(decl_asy));

cdecl_asy = cos(decl_asy);

sdecl_asy = sin(decl_asy);

crasc_asy = cos(rasc_asy);

srasc_asy = sin(rasc_asy);

% compute unit asymptote vector

shat(1) = cdecl_asy * crasc_asy;

shat(2) = cdecl_asy * srasc_asy;

shat(3) = sdecl_asy;

% compute b-plane T unit vector

that = cross(shat, pv) / norm(cross(shat, pv));

% compute b-plane R unit vector

rvec = cross(shat, that) / norm(cross(shat, that));

% sine and cosine of b-plane angle

st = sin(theta);

ct = cos(theta);

% rc = cos(theta) * rvec

rc = ct * rvec;

% ts = sin(theta) * that

ts = st * that;

% compute unit angular momentum vector
% (hhat = that * sin(theta) - rhat * cos(theta)

hhat = ts - rc;

% cosine of asymptote true anomaly

cf = -mu / (rpmag * vinf * vinf + mu);

% sine of asymptote true anomaly

sf = sqrt(1.0d0 - cf * cf);

% vtemp = hhat x shat

vtemp = cross(hhat, shat);

% vtemp = sin(ta_asy) * vtemp

vtemp = sf * vtemp;

% rphat = cosine(ta_asy) * shat

rphat = cf * shat;

% rphat = rphat - vtemp

rphat = rphat - vtemp;

% perigee velocity of launch hyperbola

vpmag = sqrt((2.0d0 * mu + rpmag * vinf * vinf) / rpmag);

% compute unit velocity vector at perigee
% (vphat = hhat x rphat)

vphat = cross(hhat, rphat) / norm(cross(hhat, rphat));

% position vector of launch hyperbola
% (rhyper = rpmag * rphat)

rhyper = rpmag * rphat;

% velocity vector of launch hyperbola
% (vhyper = vpmag * vphat)

vhyper = vpmag * vphat;

