function [fid, otype, jdtdb_tip, ddays1, jdtdb_arrival, ddays2, hp, aziml, ...
          xlatgs, rpt, xinct, fpa_tar, alt_tar, itar_type, ...
          bdott_user, bdotr_user, thetat] = e2m_readdata(filename)

% read e2m data file

% NOTE: all angular elements are returned in radians

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dtr = pi / 180.0;

% open data file

fid = fopen(filename, 'r');

% check for file open error

if (fid == -1)
    
    clc; home;
    
    fprintf('\n\n error: cannot find this file!!');
    
    pause
    
    return;
    
end

% read 70 lines of data file

for i = 1:1:70

    cline = fgetl(fid);

    switch i
        
        case 16
            
            % type of optimization

            otype = str2double(cline);
            
        case 19
            
            % launch calendar date

            tl = size(cline);

            ci = strfind(cline, ',');

            % extract month, day and year

            month = str2double(cline(1:ci(1)-1));

            day = str2double(cline(ci(1)+1:ci(2)-1));

            year = str2double(cline(ci(2)+1:tl(2)));

            jdtdb_tip = julian(month, day, year);
            
        case 22
            
            % search boundary for launch date (days)
            
            ddays1 = str2double(cline);
            
        case 25
            
            % arrival calendar date

            tl = size(cline);

            ci = strfind(cline, ',');

            % extract month, day and year

            month = str2double(cline(1:ci(1)-1));

            day = str2double(cline(ci(1)+1:ci(2)-1));

            year = str2double(cline(ci(2)+1:tl(2)));

            jdtdb_arrival = julian(month, day, year);
            
        case 28
            
            % search boundary for launch date (days)
            
            ddays2 = str2double(cline);
            
        case 35
            
            % perigee altitude of launch hyperbola

            hp = str2double(cline);
            
        case 38
            
            % launch azimuth

            aziml = dtr * str2double(cline);
            
        case 41
            
            % launch site latitude

            xlatgs = dtr * str2double(cline);
            
        case 49
            
            % type of targeting

            itar_type = str2double(cline);
            
        case 52
            
            % user-defined B dot T
            
            bdott_user = str2double(cline);
            
        case 55
            
            % user-defined B dot R
            
            bdotr_user = str2double(cline);
            
        case 58
            
            % mars periapsis radius (kilometers)

            rpt = str2double(cline);
            
        case 61
            
            % mars orbit inclination

            xinct = dtr * str2double(cline);
            
        case 64
            
            % EI flight path angle
            
            fpa_tar = dtr * str2double(cline);
            
        case 67
            
            % EI altitude
            
            alt_tar = str2double(cline);
            
        case 70
            
            % user-defined b-plane angle
            
            thetat = dtr * str2double(cline);
            
    end
    
end

% if necessary, set target flight path angle to zero

if (itar_type ~= 3)
    
    fpa_tar = 0.0;
    
end

fclose(fid);

