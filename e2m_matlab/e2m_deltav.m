function [f, g] = e2m_deltav (x)
 
% two-body, patched-conic delta-v objective function

% input

%  x = current values for launch and arrival julian dates wrt reference

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global aunit smu ip1 ip2 jdtdb0

global otype rito vito dv1 dv2

% current julian dates

jdtdb_tip = x(1) + jdtdb0;

jdtdb_arrival = x(2) + jdtdb0;

% time-of-flight

taud = jdtdb_arrival - jdtdb_tip;

tof = taud * 86400.0;

% compute initial state vector

svi = jplephem(jdtdb_tip, ip1, 11);

ri = aunit * svi(1:3);

vi = aunit * svi(4:6) / 86400.0;

% compute final state vector

svf = jplephem(jdtdb_arrival, ip2, 11);

rf = aunit * svf(1:3);

vf = aunit * svf(4:6) / 86400.0;

% solve Lambert's problem

revmax = 0;

sv1(1:3) = ri;

sv1(4:6) = vi;
    
sv2(1:3) = rf;

sv2(4:6) = vf;

[vito, ~] = glambert(smu, sv1, sv2, tof, revmax);

rito = ri;

% calculate departure delta-v

dv1(1) = vito(1) - vi(1);
dv1(2) = vito(2) - vi(2);
dv1(3) = vito(3) - vi(3);

dvm1 = norm(dv1);

% propagate transfer orbit

[~, v2] = twobody2 (smu, tof, ri, vito);

% calculate arrival delta-v

dv2(1) = vf(1) - v2(1);
dv2(2) = vf(2) - v2(2);
dv2(3) = vf(3) - v2(3);

dvm2 = norm(dv2);

% load scalar objective function

switch otype
    
case 1
    
   % launch
   
   f(1) = dvm1;
   
case 2
    
   % arrival
   
   f(1) = dvm2;
   
case 3
    
   % launch + arrival
   
   f(1) = dvm1 + dvm2;
   
case 4
    
   f(1) = dvm1 + dvm2;
   
end

% no derivatives

g = [];
