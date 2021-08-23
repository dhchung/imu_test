
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
