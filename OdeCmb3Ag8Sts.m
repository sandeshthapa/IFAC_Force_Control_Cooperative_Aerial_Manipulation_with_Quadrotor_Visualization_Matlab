close all; clear all; clc; 
tspan = [0:0.01:2.0];
% x1 = [0;-0.3;0.2];
x1 = [0.1;-0.30;0.7;0;0;-0.1;0;-pi/2];
% x2 = [0;0.3;0.2];
x2 = [0;0.30;0.68;0;0;0.1;0;pi/2];
xc = [0;0;0.5];
xcDot = [0;0;0];
McHat1 = [0.1];
McHat2 = [0.1];
q_int = [x1;x2;xc;xcDot;McHat1;McHat2];
[t,q] = ode45('Cmb3Ag8Sts',tspan,q_int);

%%

set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

fontname = 'cmss';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);

fontsize = 16;
set(0,'defaultaxesfontsize',fontsize);
set(0,'defaulttextfontsize',fontsize);
figure
plot(q(:,1),q(:,2),'-k',q(:,9),q(:,10),'-.k',q(:,17),q(:,18),':r','LineWidth',1.5)
xlabel('x position (m)')
ylabel('y postion (m)')
legend('Agent1','Agent2','Load','Location','best');
% axis([0 7 -10 1.0])
hold on 

figure
plot(t,q(:,20),'--k',t,q(:,21),':k',t,q(:,22),':r','LineWidth',1.5)
% plot(t,q(:,10))
xlabel('time (sec)')
ylabel('velocity of load (m/s)')
legend('X Velocity','Y Velocity','Z Velocity','Location','best');
hold on 

figure
plot(t,(q(:,23)),'--k',t,(q(:,24)),':k',t,0.5*ones(size(q(:,24))),':r','LineWidth',1.5)
xlabel('time (sec)')
ylabel('Mass (kg)')
% a = annotation('textarrow',[.65 .35]/2,[0.3 0.42]/0.55,'String','Actual mass of payload','Interpreter','latex');
% a.Color = 'black'; 
% a.FontSize = 16; 
legend('Estimated Mc,1','Estimated Mc,2','Actual Mc','Location','best');
axis([0 2 0 0.6])
hold on 

figure
v1 = diff(q(:,1))'./.01;
v2 = diff(q(:,2))'./.01;
v3 = diff(q(:,3))'./.01;
v9 = diff(q(:,9))'./.01;
v10 = diff(q(:,10))'./.01;
v11 = diff(q(:,11))'./.01;
    plot(tspan(1:end-1),v1,':r',tspan(1:end-1),v9,'--k',tspan(1:end-1),v2,'m',tspan(1:end-1),v10,'-.b',tspan(1:end-1),v3,':g',tspan(1:end-1),v11,'--k','LineWidth',1.5)
    xlabel('time (sec)')
    ylabel('Velocity of agents (m/s)')
    legend('Agent1,x','Agent2,x','Agent1,y','Agent2,y','Agent1,z','Agent2,z','Location','best','Orientation','vertical');
hold on

%%
figure
plot(t,q(:,1),'r',t,q(:,9),'b','LineWidth',2)
% plot(t,q(1))
xlabel('time (sec)')
ylabel('x postion (m)')
legend('Agent 1(LAM)','Agent 2(RAM)','Location','best');
hold on 

figure
plot(t,q(:,2),'r',t,q(:,10),'b','LineWidth',2)
xlabel('time (sec)')
ylabel('y postion (m)')
legend('Agent 1(LAM)','Agent 2(RAM)','Location','best');
hold on 

figure
plot(t,q(:,6),'r',t,q(:,14),'b','LineWidth',2)
xlabel('time (sec)')
ylabel('roll angle(rad)')
legend('Agent 1(LAM)','Agent 2(RAM)','Location','best');
hold on 

figure
plot(t,q(:,17),'r',t,q(:,18),'m',t,q(:,19),'b','LineWidth',2)
xlabel('time (sec)')
ylabel('position of load (m)')
legend('x pos','y pos','z pos','Location','best');
hold on 

figure
plot(t,q(:,20),'r',t,q(:,21),'m',t,q(:,22),'b','LineWidth',2)
% plot(t,q(:,10))
xlabel('time (sec)')
ylabel('velocity of load (m/s)')
legend('X velocity','Y Velocity','Z Velocity','Location','best');
hold on 

figure
plot(t,(q(:,23)+q(:,24)),'b',t,0.5*ones(size(q(:,23))),'r','LineWidth',2)
xlabel('time (sec)')
ylabel('Mass of Unknown Payload')
a = annotation('textarrow',[.65 .35]/2,[0.3 0.42]/0.55,'String','Actual mass of payload','Interpreter','latex');
a.Color = 'black'; 
a.FontSize = 16; 
hold on 

% colorspec = {[1 0 1]; [1 0 0]; [0 0 1]; ...
%   [0 1 0]};
% colorstring = 'kbgry';

figure
v1 = diff(q(:,1))'./.01;
v9 = diff(q(:,9))'./.01;
v2 = diff(q(:,2))'./.01;
v10 = diff(q(:,10))'./.01;
    plot(tspan(1:end-1),v1,'m',tspan(1:end-1),v9,'b',tspan(1:end-1),v2,'r',tspan(1:end-1),v10,'k','LineWidth',2)
    xlabel('time (sec)')
    ylabel('Velocity of agents (m/s)')
    legend('Agent 1x','Agent 2x','Agent 1y','Agent 2y','Location','best');
hold on

%% 
%Animation 
% ht = title(sprintf('Time: %0.2f sec', t(1)));
figure
n = length(q);
rotC = zeros(n,1);
xR = zeros(n,1);
yR = zeros(n,1);
zR = zeros(n,1);
pbNew1 = zeros(n,3);
rotC2 = zeros(n,1);
xR2 = zeros(n,1);
yR2 = zeros(n,1);
zR2 = zeros(n,1);
pbNew2 = zeros(n,3);
% f1 = zeros(n,1); 
% f2 = zeros(n,1);
teta1 = linspace(0,2*pi,n)';

v = VideoWriter('Videos/PayTransOct29A.mp4','MPEG-4');
open(v)
XM1Ag1  = cell(n,1);
YM1Ag1  = cell(n,1);
ZM1Ag1  = cell(n,1);

ht = title(sprintf('Time: %0.2f sec', t(1)));
for i = 1:length(q)
    
    xb1 = q(i,1);
    yb1 = q(i,2);
    zb1 = q(i,3);
    psi1  = q(i,4);      % z angle in the order ZYX [psi1 theta1 phi1] -->> [Zangle Yangle Xangle]
    theta1 = q(i,5);    %  y angle % Yaw - Pitch -Roll is not same everywhere, you can define based on your conventions
    phi1 = q(i,6);       % x angle  --->>ZYX in the order -->>YAW-PITCH-ROLL
    n11 = q(i,7);
    n21 = q(i,8);
    
    
    % Rotation matrix
    R_b1 = [ cos(psi1)*cos(theta1), cos(psi1)*sin(phi1)*sin(theta1) - cos(phi1)*sin(psi1), sin(phi1)*sin(psi1) + cos(phi1)*cos(psi1)*sin(theta1);
        cos(theta1)*sin(psi1), cos(phi1)*cos(psi1) + sin(phi1)*sin(psi1)*sin(theta1),   cos(phi1)*sin(psi1)*sin(theta1) - cos(psi1)*sin(phi1);
        -sin(theta1),         cos(theta1)*sin(phi1),                              cos(phi1)*cos(theta1)];
    %     p_b = [0 -0.3 0]';
    
    T_b1 = [ 0, -sin(psi1), cos(psi1)*cos(theta1);
        0,  cos(psi1), cos(theta1)*sin(psi1);
        1,         0,         -sin(theta1)];
    
    l1 = 0.2;
    l2 = 0.2;
   
    % Forward Kinematics
    p_b1 = [xb1 yb1 zb1]';
    pbNew1 = R_b1*p_b1;
    r= 0.15;
    %     n = length(q);
    teta = linspace(-pi,pi,n)';
    
%     x_1=r;
%     y_1=r;
%     z_1 = 0;
%     rotC = R_b1*[x_1 y_1 0]';
%     xR = rotC(1)*cos(teta);
%     yR = rotC(2)*sin(teta);
    %     zR = 0.6*ones(n,1);
    
    h1 = [1 0 0]';
    P10 = [0 0 0.0]';
    p1c = [0 0 -l1/2]';
    
    R01 = expm((hat(h1))*n11);
    
    p1c_b = P10 + R01*p1c;    % in base frame
    p1c_i = p_b1 + R_b1*p1c_b;  % in inertial frame
    
    
    % position of center of link 2
    h2 = [0 1 0]';
    p12 = [0 0 -l1]';
    p2T = [0 l2 0]';
    p2c = [0 l2/2 0]';
    
    R12 = expm((hat(h2))*n21);
    R02 = R01*R12;
    
    % postion of link 1 ;
    p1 = p_b1 + (P10 + R01*p12);
   
    % position of end effector
    p2_b1 = P10 + R01*p12 + R02*p2T; % in base frame
    x1 = p_b1 + R_b1*p2_b1;  % in inertial frame
    
    % Agent 2
    %Kinematics of right Aerial Manipulator (subscript 2)
    xb2 = q(i,9);
    yb2 = q(i,10);
    zb2 = q(i,11);
    
    psi2 = q(i,12);     % z angle (Z)
    theta2 = q(i,13);   % y angle (Y)
    phi2 = q(i,14);      % x angle (x) -->> ZYX pair yaw-pitch roll pair oder is imporant!!!
    n12 = q(i,15);
    n22 = q(i,16);
%     
    mUAV = 5.0;
    mLink1 = 0.5;
    mLink2 = 0.5;
    
    m2 = mUAV + mLink1 + mLink2;
    
    % Rotation matrix
    R_b2 = [ cos(psi2)*cos(theta2), cos(psi2)*sin(phi2)*sin(theta2) - cos(phi2)*sin(psi2), sin(phi2)*sin(psi2) + cos(phi2)*cos(psi2)*sin(theta2);
        cos(theta2)*sin(psi2), cos(phi2)*cos(psi2) + sin(phi2)*sin(psi2)*sin(theta2),   cos(phi2)*sin(psi2)*sin(theta2) - cos(psi2)*sin(phi2);
        -sin(theta2),         cos(theta2)*sin(phi2),                              cos(phi2)*cos(theta2)];
    
    % Transformation Matrix
    T_b2 = [ 0, -sin(psi2), cos(psi2)*cos(theta2);
        0,  cos(psi2), cos(theta2)*sin(psi2);
        1,         0,         -sin(theta2)];
    
    % Forward Kinematics
    p_b2 = [xb2 yb2 zb2]'; %postion of UAV with respect to inertial frame
    pbNew2 = R_b2*p_b2;
%     
%     r2= 0.15;
%     %     n = length(q);
%     teta2 = linspace(-pi,pi,n)';
%     
%     x_2=r2;
%     y_2=r2;
%     rotC2 = R_b2*[x_2 y_2 zb2]';
%     xR2 = r2*cos(teta2);
%     yR2 = r2*sin(teta2);
    %     zR = 0.6*ones(n,1);
    
    h12 = [1 0 0]';
    P102 = [0 0 0.0]';
    p1c2 = [0 0 -l1/2]';
    
    R012 = expm((hat(h12))*n12);
    
    p1c_b2 = P102 + R012*p1c2;    % in base frame
    p1c_i2 = p_b2 + R_b2*p1c_b2;  % in inertial frame
    
    % position of center of link 2
    h22 = [0 -1 0]';
    p122 = [0 0 -l1]';
    p2T2 = [0 -l2 0]';
    p2c2 = [0 -l2/2 0]';
    
    R122 = expm((hat(h22))*n22);
    R022 = R012*R122;
    
    % postion of link 1 ;
    p2 = p_b2 + (P102 + R012*p122);
    
    % position of end effector
    p2_b2 = P102 + R012*p122 + R022*p2T2; % in base frame
    
    x2 = p_b2 + R_b2*p2_b2;  % in inertial frame
    
    % Load 
    Mc = 0.5; 
    e = [0 0 1]';
    g = 9.81;
    xc = [q(i,17);q(i,18);q(i,19)];
    xcDot = [q(i,20);q(i,21);q(i,22)];
    theta = 0;
    R_theta = [cos(theta) -sin(theta) 0;
        sin(theta) cos(theta)  0;
        0          0           1];
    % xc = [0;0;0.5];
    r1 = [0 -0.10 0]';
    r2 = [0 0.10 0]';
    
    a1 = xc + R_theta*r1;
    
    a2 = xc + R_theta*r2;
    
    z1 = x1 - a1;
    % z1
    % max(z1)
    
    
    
    z2 = x2 - a2;
%     k = 100;
k = 2.5*10e4;
%     f1 = k*z1;
%     f2 = k*z2;
%     k = 100;
    f1 = k*norm(z1,3)*z1;
    f2 = k*norm(z2,3)*z2;
    
    % f_sg = 0.5*Mc*g + 10;
    f1d = [0 0.5 0.5*Mc*g]';
    f2d = [0 -0.5 0.5*Mc*g]';
    
    v1d = [.2 0 0]';
    v2d = [.2 0 0]';
%     
%     
%     v2 = [x1;p1;pbNew1];
   
    plot3(x1(1,:),x1(2,:),x1(3,:),'<',...
        p1(1,:),p1(2,:),p1(3,:),'-o',...
        p_b1(1,:),p_b1(2,:),p_b1(3,:),'o',...
        x2(1,:),x2(2,:),x2(3,:),'<',...
        p2(1,:),p2(2,:),p2(3,:),'-o',...      
        p_b2(1,:),p_b2(2,:),p_b2(3,:),'o','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','r','Color','b');
    
%         xR +p_b1(1,:),yR+p_b1(2,:),zR+ p_b1(3,:) ,...
%         xR2+p_b2(1,:),yR2 +p_b2(2,:),zR2+ p_b2(3,:) ,...
    
        line([x1(1,:),p1(1,:)],[x1(2,:),p1(2,:)],[x1(3,:),p1(3,:)],'LineWidth',3,'Color','g')
        line([p1(1,:),p_b1(1,:)],[p1(2,:),p_b1(2,:)],[p1(3,:),p_b1(3,:)],'LineWidth',3,'Color','k')
        line([x2(1,:),p2(1,:)],[x2(2,:),p2(2,:)],[x2(3,:),p2(3,:)],'LineWidth',3,'Color','g')
        line([p2(1,:),p_b2(1,:)],[p2(2,:),p_b2(2,:)],[p2(3,:),p_b2(3,:)],'LineWidth',3,'Color','k')
        
        % location of rotors
        AM1R1_b = [0.15 0 0]';
        AM1R1_i = R_b1*AM1R1_b + p_b1;
        AM2R1_i = R_b2*AM1R1_b + p_b2;
        
        AM1R2_b = [0 0.15 0]';
        AM1R2_i = R_b1*AM1R2_b + p_b1;
        AM2R2_i = R_b2*AM1R2_b + p_b2;
        
        AM1R3_b = [-0.15 0 0]';
        AM1R3_i = R_b1*AM1R3_b + p_b1;
        AM2R3_i = R_b2*AM1R3_b + p_b2;
        
        AM1R4_b = [0 -0.15 0]';
        AM1R4_i = R_b1*AM1R4_b + p_b1;
        AM2R4_i = R_b2*AM1R4_b + p_b2;
        
        % circle for motors 
%         rM1 = 0.05; 
%         xM1 = 0.05; 
%         yM1 = 0.05; 
%         
%         n = length(q);
%         
%         r1 = 0.1; 
%         x_1=0.1;
%         y_1=0.1;
%         rotC = R_b1*[x_1 y_1 0.7]';
% %         xR2 = x_1*cos(teta2);
% %         yR2 = y_2*sin(teta2);
%         zR = 0.7*ones(n,1);
%         hold on 
%         
%          rad1 = 0.1*ones(n,1); 
%         for i = i
            centerM1 = AM1R1_i;
            rad1 = 0.025;
            teta1 = linspace(0,2*pi,n)';
            xM1Ag1 = rad1*cos(teta1) + centerM1(1,:);
            yM1Ag1 = rad1*sin(teta1) + centerM1(2,:);
            zM1Ag1 = 0.00*sin(teta1)+ centerM1(3,:) + 0.015;
            %          circle1 = R_b1*[xM1Ag1';yM1Ag1;zM1Ag1];

            centerM2 = AM1R2_i;
            xM1Ag2 = rad1*cos(teta1) + centerM2(1,:);
            yM1Ag2 = rad1*sin(teta1) + centerM2(2,:);
            zM1Ag2 = 0.00*sin(teta1)+ centerM2(3,:) + 0.015;

            hold on
            cirlce1 = plot3(xM1Ag1, yM1Ag1, zM1Ag1,'LineWidth',1.5,'Color','b');
%             hold on
            plot3(xM1Ag2, yM1Ag2, zM1Ag2,'LineWidth',1.5,'Color','r')
%             hold on 
            
%         end
        
%         hold on 
        line([p_b1(1,:),AM1R1_i(1,:)],[p_b1(2,:),AM1R1_i(2,:)],[p_b1(3,:),AM1R1_i(3,:)],'LineWidth',3.5,'Color','r')
        line([p_b1(1,:),AM1R2_i(1,:)],[p_b1(2,:),AM1R2_i(2,:)],[p_b1(3,:),AM1R2_i(3,:)],'LineWidth',3.5,'Color','b')
        line([p_b1(1,:),AM1R3_i(1,:)],[p_b1(2,:),AM1R3_i(2,:)],[p_b1(3,:),AM1R3_i(3,:)],'LineWidth',3.5,'Color','r')
        line([p_b1(1,:),AM1R4_i(1,:)],[p_b1(2,:),AM1R4_i(2,:)],[p_b1(3,:),AM1R4_i(3,:)],'LineWidth',3.5,'Color','b')
        
        
        
        line([p_b2(1,:),AM2R1_i(1,:)],[p_b2(2,:),AM2R1_i(2,:)],[p_b2(3,:),AM2R1_i(3,:)],'LineWidth',3.5,'Color','r')
        line([p_b2(1,:),AM2R2_i(1,:)],[p_b2(2,:),AM2R2_i(2,:)],[p_b2(3,:),AM2R2_i(3,:)],'LineWidth',3.5,'Color','b')
        line([p_b2(1,:),AM2R3_i(1,:)],[p_b2(2,:),AM2R3_i(2,:)],[p_b2(3,:),AM2R3_i(3,:)],'LineWidth',3.5,'Color','r')
        line([p_b2(1,:),AM2R4_i(1,:)],[p_b2(2,:),AM2R4_i(2,:)],[p_b2(3,:),AM2R4_i(3,:)],'LineWidth',3.5,'Color','b')
           
    [X,Y,Z] = sphere;
    rr = 0.10;
    surface(X*rr+xc(1,:),Y*rr+xc(2,:),Z*rr+xc(3,:))
    
    axis ([-0.4 1.2 -0.4 0.4 0.2 0.8])
% axis ([-0.4 5.2 -2.4 2.4 0.2 2.0])
    xlabel(' x axis')
    zlabel('z axis')
    ylabel('y axis')
    drawnow update
    grid on
    %     axis tight manual
    M(i) = getframe(gcf);
    im = frame2im(M(1));
    %     getframe(gcf,),
    set(gca,'nextplot','replacechildren');
    title(sprintf('Time: %0.2f sec', t(i)));
%     view(az, el);
%    view([-45 145])
%     rotate3d  
%     pause(0.1)   
end 
writeVideo(v,M)
close(v)

