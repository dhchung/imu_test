function imu_3d_test()
fprintf('Loading IMU\n');
IMU = load('imu_raw.txt');
fprintf('IMU loaded\n');
fprintf('Loading GPS\n');
GPS = loadGPS('gps_raw.txt');
GPS = ProcessGPS(GPS);  %time, x, y, z, yaw

GPS(:,2:3) = GPS(:,2:3) - GPS(1,2:3);

timeline = [GPS(:,1), zeros(size(GPS,1),1), (1:size(GPS,1))';...
            IMU(:,1), ones(size(IMU,1),1), (1:size(IMU,1))'];

timeline = sortrows(timeline,1);
P = diag([0.1, 0.1, 0.1, 100, 100, 100, 100, 100]);
Q = diag([IMU_variance(1), IMU_variance(2), IMU_variance(3)]);



end