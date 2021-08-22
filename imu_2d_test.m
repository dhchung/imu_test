function imu_2d_test()
radius = 2.0;
omega = 10*pi/180.0; %rad/s

dt = 0.01;
time = 0:dt:200;
time_span = size(time,2);

inital_pose = [radius; 0; pi/2];
actual_x = radius * cos(omega * time);
actual_y = radius * sin(omega * time);
actual_psi = omega*time + pi/2;

actual_ax = 0.0;
actual_ay = radius * omega^2;
actual_w = omega;

actual_bias_x = 4;
actual_bias_y = 4;
actual_bias_w = 0.2;

actual_bias = [actual_bias_x;...
               actual_bias_y;...
               actual_bias_w];

initial_bias = [0;0;0];
% initial_uv = [omega*radius;0];
initial_uv = [0;0];

IMU_variance = [0.4, 0.4, 0.4];

GPS_variance = [0.0001, 0.0001, 0.000001];

initial_state = [inital_pose; initial_uv; initial_bias];
state_now = initial_state;
state = initial_state;

IMU_input = [repmat(actual_ax + actual_bias_x, 1, time_span);...
             repmat(actual_ay + actual_bias_y, 1, time_span);...
             repmat(actual_w + actual_bias_w, 1, time_span)];

% IMU_input = [repmat(actual_ax, 1, time_span);...
%              repmat(actual_ay, 1, time_span);...
%              repmat(actual_w, 1, time_span)];
         
         
         
GPS_input = [actual_x + randn(1, time_span)*sqrt(GPS_variance(1));...
             actual_y + randn(1, time_span)*sqrt(GPS_variance(2));...;
             actual_psi + randn(1, time_span)*sqrt(GPS_variance(3))];
         

IMU_input(1,:) = IMU_input(1,:) + randn(1, time_span)*sqrt(IMU_variance(1));
IMU_input(2,:) = IMU_input(2,:) + randn(1, time_span)*sqrt(IMU_variance(2));
IMU_input(3,:) = IMU_input(3,:) + randn(1, time_span)*sqrt(IMU_variance(3));

actual_traj = zeros(2, time_span);

P = diag([0.1, 0.1, 0.1, 100, 100, 100, 100, 100]);
Q = diag([IMU_variance(1), IMU_variance(2), IMU_variance(3)]);


H = zeros(3, size(state_now,1));
H(1:3,1:3) = eye(3);


fixed_state = [initial_state];

for i=1:time_span
    actual_pose = [actual_x(i); actual_y(i); actual_psi(i)];
    actual_traj(:,i) = actual_pose(1:2,1);
    
    figure(1);

    hold off;
    axis equal;
    
    if i>2
%         state_now = processIMU(state_now, dt, IMU_input(:,i));
        [state_now, P] = KalmanFilterIMU(P, state_now, IMU_input(:,i), Q, dt);
        state = [state state_now];
        
        if(rem(i,1)==0)
            figure(1);
            plot(actual_traj(1,1:i), actual_traj(2,1:i),'r');
            hold on;
            draw_vehicle(actual_pose, 0.2, 'b');
            plot(state(1,:), state(2,:), 'g');
            draw_vehicle(state_now(1:3,1), 0.2, 'k');
            hold off;
            axis equal;
            drawnow;
        end
        
        if(rem(i,50)==0)
            %GPS INPUT
            GPS = GPS_input(:,i);
            y = GPS - H*state_now;
            S = H*P*H' + diag(GPS_variance);
            K = P*H'*inv(S);
            state_now = state_now + K*y;
            P = (eye(size(state_now,1)) - K*H)*P;
            fixed_state = [fixed_state state_now];
            
            estimated_bias = state_now(6:8,1);
            actual_bias - estimated_bias
        end
        
    end
    
end



end

function draw_vehicle(pose, size, color)
    P1 = [1; 0]*size;
    P2 = [-1; 0.5]*size;
    P3 = [-1; -0.5]*size;

    P1 = translatePoint(P1, pose);
    P2 = translatePoint(P2, pose);
    P3 = translatePoint(P3, pose);
    
    P = [P1 P2 P3 P1];
    
    plot(P(1,:), P(2,:), color);
end

function R = rotationMat(yaw)
    R = [cos(yaw) -sin(yaw);...
         sin(yaw) cos(yaw)];
end

function T = translationMat(pose)
    x = pose(1);
    y = pose(2);
    yaw = pose(3);
    
    T = [rotationMat(yaw) [x;y];...
         0 0 1];
end

function P2 = translatePoint(P1, pose)
    T = translationMat(pose);
    Ptemp = [P1;1];
    Ptemp = T*Ptemp;
    P2 = Ptemp/Ptemp(3);
    P2 = P2(1:2,1);
end

function state2 = processIMU(state1, dt, IMU)
    xy1 = state1(1:2,1);
    psi1 = state1(3,1);
    uv1 = state1(4:5,1);
    bias1 = state1(6:8,1);
    
    state2 = state1 + [rotationMat(psi1) * uv1;...
                       IMU(3)-bias1(3);...
                       IMU(1:2,1)-bias1(1:2,1) + [uv1(2)*(IMU(3)-bias1(3));...
                                                 -uv1(1)*(IMU(3)-bias1(3))];...
                       0;...
                       0;...
                       0] * dt;
end

function A = state_transition(state, IMU, dt)
A = eye(size(state,1));
psi = state(3,1);
u = state(4,1);
v = state(5,1);
bax = state(6,1);
bay = state(7,1);
bw = state(8,1);

ax = IMU(1);
ay = IMU(2);
w = IMU(3);

A(1,3) = (-u*sin(psi)-v*cos(psi))*dt;
A(1,4) = cos(psi)*dt;
A(1,5) = -sin(psi)*dt;

A(2,3) = (u*cos(psi) - v*sin(psi))*dt;
A(2,4) = sin(psi)*dt;
A(2,5) = cos(psi)*dt;

A(3,8) = -dt;

A(4,5) = (w-bw)*dt;
A(4,6) = -dt;
A(4,8) = -v*dt;

A(5,4) = -(w-bw)*dt;
A(5,7) = -dt;
A(5,8) = u*dt;
end

function B = input_transition(state, IMU, dt)
B = zeros(size(state,1), size(IMU,1));
psi = state(3,1);
u = state(4,1);
v = state(5,1);
bax = state(6,1);
bay = state(7,1);
bw = state(8,1);

ax = IMU(1);
ay = IMU(2);
w = IMU(3);

B(3,3) = dt;
B(4,1) = dt;
B(4,3) = v*dt;
B(5,2) = dt;
B(5,3) = -u*dt;

end

function [state2, P2] = KalmanFilterIMU(P1, state1, IMU, Q, dt)
    state2 = processIMU(state1, dt, IMU);
    A = state_transition(state1, IMU, dt);
    B = input_transition(state1, IMU, dt);
    state3 = A*state1 + B*IMU;
    
    P2 = A*P1*A' + B*Q*B';
end
