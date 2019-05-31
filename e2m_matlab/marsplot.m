function marsplot()

% graphics display for e2m_matlab.m script

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global aunit jdtdb_soi jdtdb_ca rsc_soi vsc_soi

iniz = 1;

% "plot" sphere-of-influence for mars (kilometers)

rsoi_mars = 577231.618115568;

% equatorial radius of mars (kilometers)

req_mars = 3396.19;

% J2000 ecliptic-to-equatorial transformation matrix

eq2000 = [[1.000000000000000 0  0]; ...
    [0   0.917482062069182  -0.397777155931914]; ...
    [0   0.397777155931914   0.917482062069182]];

% total flight time (days)

tof = jdtdb_ca - jdtdb_soi;

fprintf('\n\nplease input the plot step size (days)\n');

dtstep = input('? ');

dtstep2 = 0.25;

fprintf('\n  please wait, computing graphics data ...\n\n');

% create heliocentric coordinate system axes vectors

hxaxisx = [0 0.1];
hxaxisy = [0 0];
hxaxisz = [0 0];

hyaxisx = [0 0];
hyaxisy = [0 0.1];
hyaxisz = [0 0];

hzaxisx = [0 0];
hzaxisy = [0 0];
hzaxisz = [0 0.1];

% create coordinate system axes vectors

xaxisx = [1 3];
xaxisy = [0 0];
xaxisz = [0 0];

yaxisx = [0 0];
yaxisy = [1 3];
yaxisz = [0 0];

zaxisx = [0 0];
zaxisy = [0 0];
zaxisz = [1 3];

% integrate the 3-body trajectory and create data points

npts = 1;

npts_mo = 0;

% load initial position and velocity vectors

for i = 1:1:3
    
    yi(i) = rsc_soi(i);

    yi(i + 3) = vsc_soi(i);
    
end

% number of differential equations

neq = 6;

% rkf78 tolerance

tetol = 1.0e-10;

% create initial data points

rtmp = eq2000' * rsc_soi';

x2(npts) = rtmp(1);

y2(npts) = rtmp(2);

z2(npts) = rtmp(3);

svmars = jplephem(jdtdb_soi, 4, 11);

rmars = svmars(1:3);

rtmp = eq2000' * rmars;

x3(1) = rtmp(1);

y3(1) = rtmp(2);

z3(1) = rtmp(3);

svearth = jplephem(jdtdb_soi, 3, 11);

rearth = svearth(1:3);

rtmp = eq2000' * rearth;

x4(1) = rtmp(1);

y4(1) = rtmp(2);

z4(1) = rtmp(3);

xtime(npts) = 0.0;

ti = 0.0;

while(1)
    
    % step size guess (days)

    h = 300.0 / 86400.0;

    % increment initial and final times

    tf = ti + dtstep;

    % check for last step size

    if (tf > tof)
        
        tf = tof;
        
    end

    % integrate from ti to tf

    yfinal = rkf78('e2m_eqm3', neq, ti, tf, h, tetol, yi);

    % create transfer orbit graphics data (x2, y2, z2)

    npts = npts + 1;

    % reset initial time and state

    ti = tf;

    yi = yfinal;

    xtime(npts) = tf;

    rtmp = eq2000' * yfinal(1:3)';

    x2(npts) = rtmp(1);

    y2(npts) = rtmp(2);

    z2(npts) = rtmp(3);

    svmars = jplephem(jdtdb_soi + tf, 4, 11);

    rmars = svmars(1:3);

    % create data for Mars orbit (x3, y3, z3)

    rtmp = eq2000' * rmars;

    x3(npts) = rtmp(1);

    y3(npts) = rtmp(2);

    z3(npts) = rtmp(3);

    % create data for Earth orbit (x4, y4, z4)

    svearth = jplephem(jdtdb_soi + tf, 3, 11);

    rearth = svearth(1:3);

    rtmp = eq2000' * rearth;

    x4(npts) = rtmp(1);

    y4(npts) = rtmp(2);

    z4(npts) = rtmp(3);

    rm2sc(1) = aunit * (yfinal(1) - rmars(1));

    rm2sc(2) = aunit * (yfinal(2) - rmars(2));

    rm2sc(3) = aunit * (yfinal(3) - rmars(3));

    if (norm(rm2sc) <= rsoi_mars && iniz == 1)
        
        % reduce plot step size within Mars soi

        if (iniz == 1)
            
            dtstep = dtstep2 / 24.0;
            
            iniz = 0;
            
        end
        
    end
    
    if (norm(rm2sc) <= 50000.0)
        
        % create trajectory within 50,000 kilometers of mars (x5, y5, z5)

        jdate = jdtdb_soi + tf;

        tmatrix = mme2000(jdate);

        rtmp = tmatrix * rm2sc';

        npts_mo = npts_mo + 1;

        x5(npts_mo) = rtmp(1) / req_mars;

        y5(npts_mo) = rtmp(2) / req_mars;

        z5(npts_mo) = rtmp(3) / req_mars;
        
    end

    % check for end of simulation

    if (tf >= tof)
        
        break;
        
    end
    
end

% plot transfer trajectory

figure(1);

hold on;

grid on;

% plot transfer trajectory

plot3(x2(1), y2(1), z2(1), '.r');

plot3(x2, y2, z2, '-r');

plot3(x2(end), y2(end), z2(end), '.r');

plot3(x3, y3, z3, '-g');

plot3(x4, y4, z4, '-b');

% plot sun and coordinate system axes

plot3(0, 0, 0, '*y');

plot3(hxaxisx, hxaxisy, hxaxisz, '-r');

plot3(hyaxisx, hyaxisy, hyaxisz, '-g');

plot3(hzaxisx, hzaxisy, hzaxisz, '-b');

% label plot and axes

xlabel('X coordinate (AU)', 'FontSize', 12);

ylabel('Y coordinate (AU)', 'FontSize', 12);

zlabel('Z coordinate (AU)', 'FontSize', 12);

title('Heliocentric Transfer Trajectory', 'FontSize', 16);

axis equal;

view(-50, 20);

print -depsc -tiff -r300 marsplot1.eps;

% open second figure for areocentric view

figure(2);

hold on;

grid on;

rotate3d on;

% create martian surface

[x, y, z] = sphere(12);

h = surf(x, y, z);

colormap([0.5 0.5 0.5]);

set (h, 'edgecolor', [1 1 1]);

% plot coordinate system axes

plot3(xaxisx, xaxisy, xaxisz, '-r');

plot3(yaxisx, yaxisy, yaxisz, '-g');

plot3(zaxisx, zaxisy, zaxisz, '-b');

% plot initial position

plot3(x5(1), y5(1), z5(1), '*b');

% plot trajectory

plot3(x5, y5, z5, '-b');

% plot final position

plot3(x5(end), y5(end), z5(end), '.r');

xlabel('X coordinate (MR)', 'FontSize', 12);

ylabel('Y coordinate (MR)', 'FontSize', 12);

zlabel('Z coordinate (MR)', 'FontSize', 12);

title('Mars-centered Trajectory', 'FontSize', 16);

axis equal;

view(-50, 20);

print -depsc -tiff -r300 marsplot2.eps;


