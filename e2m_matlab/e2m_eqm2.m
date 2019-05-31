function ydot = e2m_eqm2 (t, y)

% first-order heliocentric equations of motion

% point-mass gravity of the sun, Mercury, Venus, Earth, Mars,
% Jupiter, Saturn and Uranus

% Battin's f(q) formulation

% input

%  t = time since soi (days)
%  y = spacecraft heliocentric state vector (au & au/day)

% output

%  ydot = integration vector

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rp = zeros(7, 3);

rp2sc = zeros(7, 3);

accp = zeros(3, 1);

q = zeros(3, 1);

f = zeros(3, 1);

d3 = zeros(3, 1);

global xmu jdtdb_soi

% current julian date (relative to soi event)

jdate = jdtdb_soi + t;

% distance from sun to spacecraft (au)

rsun2sc = norm(y(1:3));

rrsun2sc = -xmu(1) / rsun2sc^3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate planetary position vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1:7

    svplanet = jplephem(jdate, i, 11);

    rplanet = svplanet(1:3);

    for j = 1:1:3
        
        rp(i, j) = rplanet(j);
        
    end
    
end

% compute planet-centered position vectors of spacecraft (au)

for i = 1:1:7
    
    rp2sc(i, 1) = y(1) - rp(i, 1);
    
    rp2sc(i, 2) = y(2) - rp(i, 2);
    
    rp2sc(i, 3) = y(3) - rp(i, 3);
    
end

% compute f(q) functions for each planet

for k = 1:1:7
    
    q(k) = dot(y(1:3), y(1:3) - 2.0 * rp(k, :)') / dot(rp(k, :), rp(k, :));

    f(k) = q(k) * ((3.0 + 3.0 * q(k) + q(k) * q(k)) / (1.0 + (1.0 + q(k))^1.5));

    d3(k) = norm(rp2sc(k, :)) * norm(rp2sc(k, :)) * norm(rp2sc(k, :));
    
end

% compute planetary perturbations

for j = 1:1:3
    
    accp(j) = 0.0;

    for k = 1:1:7
        
        accp(j) = accp(j) - xmu(k + 1) * (y(j) + f(k) * rp(k, j)) / d3(k);
        
    end
    
end

% compute integration vector

ydot = [ y(4)
    y(5)
    y(6)
    accp(1) + y(1) * rrsun2sc
    accp(2) + y(2) * rrsun2sc
    accp(3) + y(3) * rrsun2sc];

