function [value, isterminal, direction] = e2m_fpa_event(t, y)

% flight path angle event function

% input

%  t = time since soi (days)
%  y = spacecraft heliocentric state vector (au, au/day)

% output

%  value = mars-centered flight path angle (radians)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jdtdb_soi fpa_tar

% mars heliocentric state vector at current time t

jdate = jdtdb_soi + t;

svmars = jplephem(jdate, 4, 11);

rmars = svmars(1:3);

vmars = svmars(4:6);

% form the mars-centered spacecraft position and velocity vectors

rm2sc = y(1:3) - rmars(1:3);

vm2sc = y(4:6) - vmars(1:3);

tmatrix = mme2000(jdate);

rsc = tmatrix * rm2sc;

vsc = tmatrix * vm2sc;

% flight path angle (radians)

fpa = dot(rsc', vsc) / (norm(rsc) * norm(vsc));

value = fpa - fpa_tar;

isterminal = 1;

direction =  [];

