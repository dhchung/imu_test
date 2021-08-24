function imu_gps_test_nobias()

fprintf('Loading IMU\n');
% IMU = load('imu_raw.txt');

IMU = load('imu.txt');
% Quaternion xyzw to wxyz
IMU(:,2:5) = [IMU(:,5) IMU(:,2:4)];

fprintf('IMU loaded\n');
fprintf('Loading GPS\n'); 
GPS = loadGPS('gps_raw.txt');
GPS = ProcessGPS(GPS);  %time, x, y, z, yaw
GPS_raw = GPS;
GPS(:,5) = GPS(:,5)/180.0*pi;


GPS(:,2:3) = GPS(:,2:3) - GPS(1,2:3);

timeline = [GPS(:,1), zeros(size(GPS,1),1), (1:size(GPS,1))';...
            IMU(:,1), ones(size(IMU,1),1), (1:size(IMU,1))'];

timeline = sortrows(timeline,1);

BX = [-1.4264;-1.0078;0.8167];

BX = [0;0;0];
BW = [0;0;0];
BA = [0;0;0];

initial_state = [0;0;0;... %x y z
                 0;0;0;... %roll pitch yaw
                 0;0;0;... %u v w
                 BX;
                 BW;
                 BA];

Bias_uncertaintyX = [10, 10, 10];
Bias_uncertaintyW = [0.4, 0.4, 0.4];
Bias_uncertaintyA = [0.4, 0.4, 0.4];

IMU_uncertaintyX = [10, 10, 10];
IMU_uncertaintyW = [1, 1, 1];
IMU_uncertaintyA = [0.1, 0.1, 0.1];

Q = diag([IMU_uncertaintyX, IMU_uncertaintyW]);

GPS_uncertaintyX = [0.5, 0.5, 0.5];
GPS_uncertaintyA = 0.1;

P = diag([0.1, 0.1, 0.1,...
          0.1, 0.1, 0.1,...
          2, 2, 2,...
          Bias_uncertaintyX,...
          Bias_uncertaintyW,...
          Bias_uncertaintyA]);

start = false;

state =[];
state_now = initial_state;

last_time = 0;

state_gps = [];

for i=1:size(timeline,1)
    if and(~start, timeline(i,2)==1)
        continue;
    end
    if and(~start, timeline(i,2)==0)
        gps_idx = timeline(i,3);
        
        start = true;
        initial_state = [GPS(gps_idx,2:4)'; [0;0;GPS(gps_idx,5)]; zeros(3,1);BX;BW;BA];
        GPS_input = GPS(gps_idx,2:5)';
        state_gps = [state_gps GPS_input];
        
        state_now = initial_state;
        state = [state state_now];
        last_time = timeline(i,1);

        continue;
    end
    
    imu_processed = false;
    
    if timeline(i,2)==1
        %process imu
        imu_idx = timeline(i,3);
        IMU_data = [IMU(imu_idx, 27:29)';...
                    IMU(imu_idx, 15:17)'];
        
        dt = timeline(i,1) - last_time;
        last_time = timeline(i,1);
        
        [state_now, P] = IMUPrediction(state_now, P, IMU_data, Q, dt);
        
        
        IMU_attitude = quat2eul(IMU(imu_idx, 2:5))';
%         [state_now, P] = IMUUpdate(state_now, P, IMU_attitude, diag(IMU_uncertaintyA));        
        
        state = [state state_now];
%         P(1:3,1:3)
        state_now(10:18)

        imu_processed = true;
        
        continue;
    end

    if timeline(i,2)==0
        %process GPS
        gps_idx = timeline(i,3);
        if(GPS_raw(gps_idx,5)==180)
            continue;
        end
            
        
        dt = timeline(i,1) - last_time;
        last_time = timeline(i,1);
        [state_now, P] = IMUPrediction(state_now, P, zeros(6,1), zeros(6), dt);
        
        GPS_input = GPS(gps_idx,2:5)';
        state_gps = [state_gps GPS_input];

        [state_now, P] = GPSUpdate(state_now, P, GPS_input, diag([GPS_uncertaintyX, GPS_uncertaintyA]));
        state = [state state_now];
        
        figure(2);
        plot3(state(1,:), state(2,:), state(3,:),'g');
        hold on;
        plot3(state_gps(1,:), state_gps(2,:), state_gps(3,:),'m');
        draw_vehicle(state_now);
        hold off;
        axis equal;
        xlabel X;
        ylabel Y;
        zlabel Z;
        view(0,90);
        drawnow;

        continue;
    end    

end

figure(1);
plot3(state(1,:), state(2,:), state(3,:));
axis equal;

plotgps(GPS);

end

function [state2, P2] = IMUPrediction(state1, P1, IMU_input, Q, dt)
    state2 = Prediction(state1, IMU_input, dt);
    hs = 0.01*ones(18,1);
    hi = 0.01*ones(6,1);
    [A, B] = NumJacobi(state1, IMU_input, hs, hi, dt);
    
    state2(4:6,1) = pi2pi(state2(4:6,1));
    
    P2 = A*P1*A' + B*Q*B';
end

function state2 = Prediction(state1, IMU_input, dt)
    Gravity = 9.8;

    IMU_A = IMU_input(1:3,1);
    IMU_W = IMU_input(4:6,1);

    RPY = state1(4:6,1);
    UVW = state1(7:9,1);
    BA = state1(10:12,1);
    BW = state1(13:15,1);

    state2 = state1 + [RotationR(RPY)*UVW;...
                       JacobianR(RPY)*(IMU_W - BW);...
                       (IMU_A-BA) + RotationR(RPY)'*[0;0;Gravity] + cross(UVW, [IMU_W-BW]);...
                       zeros(9,1)]*dt;
end

function [state2, P2] = IMUUpdate(state1, P1, IMU_attitude, IMU_var)
    H = [zeros(3) eye(3) zeros(3) zeros(3) zeros(3) -eye(3)];
    y = IMU_attitude - H*state1;
    y = pi2pi(y);
    S = H*P1*H' + IMU_var;
    K = P1*H'*inv(S);
    
    state2 = state1 + K*y;
    state2(4:6,1) = pi2pi(state2(4:6,1));
    P2 = (eye(size(state1,1)) - K*H)*P1;
end


function [state2, P2] = GPSUpdate(state1, P1, GPS_input, GPS_var)
    H = zeros(4,size(state1,1));
    H(1:3,1:3) = eye(3);
    H(4,6) = 1;
    
    y = GPS_input - H*state1;
    y(4) = pi2pi(y(4));
    S = H*P1*H' + GPS_var;
    K = P1*H'*inv(S);
    
    state2 = state1 + K*y;
    state2(4:6,1) = pi2pi(state2(4:6,1));
    P2 = (eye(size(state1,1)) - K*H)*P1;
end



function [A, B] = NumJacobi(state, imu, hs, hi, dt)
    A = zeros(size(state,1));
    B = zeros(size(state,1), size(imu,1));
    for i=1:size(A,1)
        for j=1:size(A,2)
            s_1 = state;
            s_1(j) = s_1(j) + hs(j);
            s_2 = state;
            s_2(j) = s_2(j) - hs(j);
            s_p = Prediction(s_1, imu, dt);
            s_m = Prediction(s_2, imu, dt);
            A(i,j) = (s_p(i,1)-s_m(i,1))/(2*hs(j));
        end
    end

    for i=1:size(B,1)
        for j=1:size(B,2)
            i_1 = imu;
            i_1(j) = i_1(j)+hi(j);
            i_2 = imu;
            i_2(j) = i_2(j)-hi(j);

            s_p = Prediction(state, i_1, dt);
            s_m = Prediction(state, i_2, dt);

            B(i,j) = (s_p(i,1)-s_m(i,1))/(2*hi(j));
        end
    end
end

function plotgps(GPS)
    figure(2);
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

function angle2 = pi2pi(angle)
    i = floor(sign(angle) .* angle/(2*pi));
    angle = angle - 2*pi*i.*sign(angle);
    angle(abs(angle)>pi) = angle(abs(angle)>pi) - sign(angle(abs(angle)>pi))*2*pi;
    angle2 = angle;
end


function draw_vehicle(state)
XYZ = state(1:3,1);
RPY = state(4:6,1);
P0 = XYZ;
R = RotationR(RPY);
P1 = [1; 0; 0];
P2 = [0; 1; 0];
P3 = [0; 0; 1];

P1 = R*P1;
P2 = R*P2;
P3 = R*P3;

P1 = P0+P1;
P2 = P0+P2;
P3 = P0+P3;


PX = [P0 P1];
PY = [P0 P2];
PZ = [P0 P3];

plot3(PX(1,:), PX(2,:), PX(3,:),'r');
plot3(PY(1,:), PY(2,:), PY(3,:),'g');
plot3(PZ(1,:), PZ(2,:), PZ(3,:),'b');

end