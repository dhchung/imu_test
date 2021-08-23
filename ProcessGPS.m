function GPS_new = ProcessGPS(GPS)
    UTMXY = ToUTM(GPS(2:3,:));
    GPS_new = [GPS(1,:); UTMXY; GPS(4:5,:)];
    GPS_new = GPS_new';
end