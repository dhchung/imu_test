function imu_test()
fprintf("LOADING\n");
A = load('/mnt/sdb/Avikus/20210701/exp_night_1/imu/imu.txt');
% A = load('/home/dongha/imu_raw.txt');
fprintf("LOADED\n");

Orientation = A(:,2:5);
AngularVelocity = A(:,15:17);
LinearAcceleration = A(:,27:29);
Time = A(:,1);

Gravity = 9.80655;

rpy = [0;0;0];
xyz = [0;0;0];
uvw = [0;0;0];

X = [];

for i=2:size(A,1)
    dt = Time(i,1) - Time(i-1,1);
    AV = AngularVelocity(i,:);
    LA = LinearAcceleration(i,:);
    
    dxyzdot = RotationR(rpy)*uvw;
    drpydot = JacobianR(rpy)*AV';
%     duvwdot = LA' + RotationR(rpy)*[0;0;Gravity] + cross(uvw, AV');
    duvwdot = LA' + cross(uvw, AV');

    dxyz = dxyzdot * dt;
    drpy = drpydot * dt;
    duvw = duvwdot * dt;
    
    xyz = xyz + dxyz;
    rpy = rpy + drpy;
    uvw = uvw + duvw;
    
    X = [X, [xyz;rpy;uvw]];
    
    T = [RotationR(rpy) xyz;0 0 0 1];
    Temp = [0; 100; 0; 1];
    Temp = T*Temp;
    
    Arrow = zeros(3,2);   
    Arrow(:,1) = xyz;
    Arrow(:,2) = Temp(1:3,1);
    
    
    if(rem(i, 100)==0)
        figure(1);
        plot3(X(1,:), X(2,:), X(3,:), 'b');
        hold on;
        plot3(Arrow(1,:), Arrow(2,:), Arrow(3,:), 'r');
        hold off;
        axis equal;
        xlabel X;
        ylabel Y;
        zlabel Z;
        
        lims = [xyz-10 xyz+10];
        
        axis([lims(1,:), lims(2,:), lims(3,:)]);
        drawnow;
        
    end
    
end

end


function M = JacobianR(rpy)
    roll = rpy(1,1);
    pitch = rpy(2,1);
    yaw = rpy(3,1);
    M = [1, sin(roll) * tan(pitch), cos(roll)*tan(pitch);...
         0, cos(roll), -sin(roll);...
         0, sin(roll)/cos(pitch), cos(roll)/cos(pitch)];
end

function R = RotationR(rpy)
    roll = rpy(1,1);
    pitch = rpy(2,1);
    yaw = rpy(3,1);
    Rx = [1, 0, 0;...
          0, cos(roll), -sin(roll);...
          0, sin(roll), cos(roll)];
    Ry = [cos(pitch), 0, sin(pitch);...
          0, 1, 0;...
          -sin(pitch), 0, cos(pitch)];
    Rz = [cos(yaw), -sin(yaw), 0;...
          sin(yaw), cos(yaw), 0;...
          0, 0, 1];

    R = Rz*Ry*Rx;
end