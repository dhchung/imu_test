function imu_test()
fprintf("LOADING\n");
% A = load('imu.txt');
A = load('imu_raw.txt');
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
    duvwdot = LA' + RotationR(rpy)'*[0;0;Gravity] + cross(uvw, AV');
%     duvwdot = LA' + cross(uvw, AV');

    dxyz = dxyzdot * dt;
    drpy = drpydot * dt;
    duvw = duvwdot * dt;
    
    xyz = xyz + dxyz;
    rpy = rpy + drpy;
    uvw = uvw + duvw;
    
    X = [X, [xyz;rpy;uvw]];
    
    T = [RotationR(rpy) xyz;0 0 0 1];
    Temp = [10; 0; 0; 1];
    Temp = T*Temp;
    
    Arrow = zeros(3,2);   
    Arrow(:,1) = xyz;
    Arrow(:,2) = Temp(1:3,1);
    
    
    if(rem(i, 100)==1)
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
        
%         axis([lims(1,:), lims(2,:), lims(3,:)]);
        drawnow;
        
    end
    
end

end