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
