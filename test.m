function test()

fprintf('Loading IMU\n');
% IMU = load('imu_raw.txt');
IMU = load('imu.txt');

fprintf('IMU loaded\n');
fprintf('Loading GPS\n');
GPS = loadGPS('gps_raw.txt');
GPS = ProcessGPS(GPS);  %time, x, y, z, yaw
GPS(:,5) = GPS(:,5)/180.0*pi;
GPS(:,2:3) = GPS(:,2:3) - GPS(1,2:3);
timeline = [GPS(:,1), zeros(size(GPS,1),1), (1:size(GPS,1))';...
            IMU(:,1), ones(size(IMU,1),1), (1:size(IMU,1))'];
timeline = sortrows(timeline,1);

GPS(:,5) = pi2pi(GPS(:,5));

QWXYZ = [IMU(:,5), IMU(:,2:4)];

imu_ang = quat2eul(QWXYZ);

imu_ang = [pi2pi(imu_ang(:,1)), pi2pi(imu_ang(:,2)), pi2pi(imu_ang(:,3))];
% 
% for i = 1 :size(IMU,1)
%     imu_R = quat2rotm(QWXYZ(i,:));
%     P0 = [0; 0; 0];
%     P1 = [1; 0; 0];
%     P2 = [0; 1; 0];
%     P3 = [0; 0; 1];
%     
%     P1 = imu_R*P1;
%     P2 = imu_R*P2;
%     P3 = imu_R*P3;
%     
%     PX = [P0 P1];
%     PY = [P0 P2];
%     PZ = [P0 P3];
%     
%     figure(1);
%     plot3(PX(1,:), PX(2,:), PX(3,:),'r');
%     hold on;
%     plot3(PY(1,:), PY(2,:), PY(3,:),'g');
%     plot3(PZ(1,:), PZ(2,:), PZ(3,:),'b');
%     hold off;
%     
%     axis equal;
%     xlabel X;
%     ylabel Y;
%     zlabel Z;
%     grid on;
%     
%     
%     
% end
% 



% for i = 1:size(GPS,1)
%     XY = GPS(1:i,2:3);
%     yaw = GPS(i,5);
%     
%     P1 = XY(i,1:2);
%     P2 = XY(i,1:2)' + getR(yaw)*[10;0];
%     P = [P1;P2'];
%     
%     figure(2);
%     plot(XY(:,1), XY(:,2));
%     hold on;
%     plot(P(:,1),P(:,2),'r');
%     hold off;
%     axis equal;
%     xlabel X;
%     ylabel Y;
%     drawnow;
%     
% end


GPS_LOG = [GPS(:,1), GPS(:,5)*180/pi];
IMU_LOG = [IMU(:,1), imu_ang*180/pi];

figure(1);
plot(GPS_LOG(:,1), GPS_LOG(:,2));
hold on;

plot(IMU_LOG(:,1), IMU_LOG(:,2));
plot(IMU_LOG(:,1), IMU_LOG(:,3));
plot(IMU_LOG(:,1), IMU_LOG(:,4));
hold off;

legend('GPS', 'IMU yaw', 'IMU pitch', 'IMU roll');

end

function angle2 = pi2pi(angle)
    i = floor(sign(angle) .* angle/(2*pi));
    angle = angle - 2*pi*i.*sign(angle);
    angle(abs(angle)>pi) = angle(abs(angle)>pi) - sign(angle(abs(angle)>pi))*2*pi;
    angle2 = angle;
end

function R = getR(yaw)
    R = [cos(yaw), -sin(yaw);...
         sin(yaw), cos(yaw)];

end