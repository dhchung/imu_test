function imu_gps_test()

fprintf('Loading IMU\n');
IMU = load('/mnt/sdb/Avikus/20210701/exp_night_1/imu/imu.txt');
fprintf('IMU loaded\n');
fprintf('Loading GPS\n');
GPS = loadGPS('/mnt/sdb/Avikus/20210701/exp_night_1/gps/gps_raw.txt');
GPS = ProcessGPS(GPS);  %time, x, y, z, yaw

timeline = [GPS(:,1), zeros(size(GPS,1),1), (1:size(GPS,1))';...
            IMU(:,1), ones(size(IMU,1),1), (1:size(IMU,1))'];

timeline = sortrows(timeline,1);

initial_state = [0;0;0;... %x       y       z
                 0;0;0;... %roll    pitch   yaw
                 0;0;0;... %u       v       w
                 0;0;0;... %bx      by      bz
                 0;0;0];   %br      bp      by

start = false;

for i=1:size(timeline,1)
    if and(~start, timeline(i,2)==1)
        continue;
    end
    if and(~start, timeline(i,2)==0)
        gps_idx = timeline(i,3);
        
        start = true;
        initial_state = [GPS(gps_idx,2:4)'; [0;0;GPS(gps_idx,5)]; zeros(6,1)];
        continue;
    end

    if timeline(i,2)==0
        %process GPS
        
        continue;
    end
    
    
    if timeline(i,2)==1
        %process imu
        continue;
    end
    

end
        
        
        

plotgps(GPS);

end

function plotgps(GPS)
    figure(1);
    plot3(GPS(:,2), GPS(:,3), GPS(:,4),'b');
    hold on;
    scatter3(GPS(1,2), GPS(1,3), GPS(1,4), 20, 'r');
    scatter3(GPS(end,2), GPS(end,3), GPS(end,4), 20, 'b');
    hold off;
    axis equal;
    zlim([-10 10]);
    xlabel('North X[m]');
    ylabel('East Y[m]');
    zlabel('Down Z[m]');
    set(gca, 'YDir','reverse')
    set(gca, 'ZDir','reverse')
    
end

function GPS = loadGPS(filedir)

b_GNRMC = false;
b_GNGGA = false;

GPS = [];
latitude = 0;
longitude = 0;
altitude = 0;
yaw = 0;
fix = 0;
HDOP = 0;
time = 0;

fid = fopen(filedir);
tline = fgetl(fid);
while(ischar(tline))
    C = strsplit(tline);
    datatype = C{2};
    if(strcmp(datatype, '$GNRMC'))
        b_GNRMC = true;
        b_GNGGA = false;
        

    end
    
    if(strcmp(datatype, '$GNGGA'))
        if(~b_GNRMC)
            b_GNGGA = false;
            tline = fgetl(fid);
            continue;
        end
        b_GNGGA = true;
        latitude = str2double(C{4});
        longitude = str2double(C{5});
        altitude = str2double(C{8});
        fix = str2double(C{6});
        HDOP = str2double(C{7});
    
    end
    
    if(strcmp(datatype, '$GNHDT'))
       if(~b_GNGGA)
           b_GNRMC = false;
           tline = fgetl(fid);
           continue;
       end
       if(~b_GNRMC)
           b_GNGGA = false;
           tline = fgetl(fid);
           continue;
       end
       time = str2double(C{1});
       yaw = str2double(C{3});
       GPS = [GPS [time;latitude;longitude;altitude;yaw;fix;HDOP]];
    end
    tline = fgetl(fid);
end



end


function GPS_new = ProcessGPS(GPS)

UTMXY = ToUTM(GPS(2:3,:));

GPS_new = [GPS(1,:); UTMXY; GPS(4:5,:)];

GPS_new = GPS_new';
end


function UTMXY = ToUTM(LATLON)
    UTM_ZON = 52;
    Wa = 6378137;
    Wb = 6356752.314245;
    We = 0.081819190842965;
    Weps = 0.006739496742333;

    dlat = LATLON(1,:);
    dlon = LATLON(2,:);
    
	% coordinates in radians
    lat = dlat*pi/180.0;
    lon = dlon*pi/180.0;
    % UTM parameters
    lon0_f = floor(dlon/6)*6+3; %reference longitude in degrees
    lon0 = lon0_f * pi/180.0;   %in radians
    k0 = 0.9996;                %scale on central meridian
    
    FE = 500000;				% false easting
	FN = (dlat < 0) * 10000000; % false northing

	% Equations parameters
	% N: radius of curvature of the earth perpendicular to meridian plane
	% Also, distance from point to polar axis
	WN = Wa ./ sqrt(1 - pow(We, 2) * pow(sin(lat), 2));
	WT = pow(tan(lat), 2);
	WC = (pow(We, 2) ./ (1 - pow(We, 2))) * pow(cos(lat), 2);
	WLA = (lon - lon0) .* cos(lat);
	% M: true distance along the central meridian from the equator to lat
	WM = Wa * ((1 - pow(We,2)/4 - 3*pow(We,4)/64 - 5*pow(We,6)/256)*lat...
                -(3*pow(We,2)/8 + 3*pow(We,4)/32 + 45*pow(We,6)/1024)*sin(2*lat)...
                +(15*pow(We,4)/256 + 45*pow(We,6)/1024)*sin(4*lat)...
                -(35*pow(We,6)/3072)*sin(6*lat));

	% northing
	% M(lat0) = 0 so not used in following formula
	rutmx = FN +...
            k0*WM +...
            k0*WN.*tan(lat).*(pow(WLA,2)/2 +...
            (5 - WT + 9*WC + 4*pow(WC,2)).*pow(WLA, 4)/24 +...
            (61 - 58*WT + pow(WT,2) + 600*WC - 330*Weps).*pow(WLA,6)/720);

	% easting
	rutmy = FE + ...
            k0*WN.*(WLA + (1 - WT + WC).*pow(WLA,3)/6 + (5 - 18*WT + pow(WT,2) + 72*WC - 58*Weps).*pow(WLA,5)/120);
    
    UTMXY = [rutmx;rutmy];

end

function pow_result = pow(a, b)
pow_result = a.^b;
end