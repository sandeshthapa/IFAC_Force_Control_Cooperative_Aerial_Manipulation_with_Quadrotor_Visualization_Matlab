function zDot = Cmb3Ag8Sts(t,q)

%Kinematics of Left Aerial Manipulator
xb1 = q(1); 
yb1 = q(2);
zb1 = q(3);
psi1 = q(4);      % z angle in the order ZYX [psi1 theta1 phi1] -->> [Zangle Yangle Xangle]
                  % Yaw - Pitch -Roll is not same everywhere, you can defined
                  % based on your conventions
theta1 = q(5);    % y angle
phi1 = q(6);       % x angle  --->>ZYX in the order -->>YAW-PITCH-ROLL
n11 = q(7); 
n21 = q(8);

% x_dot1 = 0.2;
% y_dot1 = 0.2;
% z_dot1 = 0;
% phi_dot1 = -0.4;  % x angle 
% theta_dot1 = 0;   % y anlge 
% psi_dot1 = 0;     % z angle  in the order ZYX  [psiDot thetaDot phiDot] -->> [z y x]Dot
% n1_dot1 = 0;
% n2_dot1 = 0;

mUAV = 5.0; 
mLink1 = 0.5; 
mLink2 = 0.5;
l1 = 0.2; 
l2 = 0.2;
 
m1 = mUAV + mLink1 + mLink2; 

% Rotation matrix
R_b1 = [ cos(psi1)*cos(theta1), cos(psi1)*sin(phi1)*sin(theta1) - cos(phi1)*sin(psi1), sin(phi1)*sin(psi1) + cos(phi1)*cos(psi1)*sin(theta1); 
      cos(theta1)*sin(psi1), cos(phi1)*cos(psi1) + sin(phi1)*sin(psi1)*sin(theta1),   cos(phi1)*sin(psi1)*sin(theta1) - cos(psi1)*sin(phi1);
      -sin(theta1),         cos(theta1)*sin(phi1),                              cos(phi1)*cos(theta1)];

  
T_b1 = [ 0, -sin(psi1), cos(psi1)*cos(theta1); 
        0,  cos(psi1), cos(theta1)*sin(psi1); 
        1,         0,         -sin(theta1)];

% Forward Kinematics
p_b = [xb1 yb1 zb1]'; %postion of UAV with respect to inertial frame
% p_b_dot1 = [x_dot1 y_dot1 z_dot1]';
% n_dot1 = [n1_dot1 n2_dot1]'; % joint angles rate 
h1 = [1 0 0]'; 
P10 = [0 0 0.0]';
p1c = [0 0 -l1/2]';

R01 = expm((hat(h1))*n11);

p1c_b = P10 + R01*p1c;    % in base frame 
p1c_i = p_b + R_b1*p1c_b;  % in inertial frame

% position of center of link 2 
h2 = [0 1 0]';
p12 = [0 0 -l1]';
p2T = [0 l2 0]';
p2c = [0 l2/2 0]';

R12 = expm((hat(h2))*n21);
R02 = R01*R12;

% position of end effector 
p2_b1 = P10 + R01*p12 + R02*p2T; % in base frame
% disp(p2_b1);
% disp(x1);
% disp(y1);
% disp(z1);
x1 = p_b + R_b1*p2_b1;  % in inertial frame
% disp('x1'); 
% disp(x1);

% J2c = [h1 h2; hat(h1)*(R01*p12) hat(h2)*p2c_b];
% Jacobian for end effector 
JT1 = [h1 h2; hat(h1)*R01*p12 hat(h2)*p2_b1];
JwT1 = JT1(1:3,1:2);
JvT1 = JT1(4:6,1:2);



% w_b1 = [phi_dot1*cos(psi1)*cos(theta1) - theta_dot1*sin(psi1); 
%       theta_dot1*cos(psi1) + phi_dot1*cos(theta1)*sin(psi1);
%                      psi_dot1 - phi_dot1*sin(theta1)];
% x1Dot = p_b_dot1 + hat(w_b1)*R_b1*p2_b1 + R_b1*JvT1*n_dot1;    % end effector


%% Agent 2 
% m2 = 1.0;
% x2 = [q(4);q(5);q(6)];
% x2Dot = [q(13);q(14);q(15)];
%Kinematics of right Aerial Manipulator (subscript 2)
xb2 = q(9); 
yb2 = q(10);
zb2 = q(11);

psi2 = q(12);     % z angle (Z)
theta2 = q(13);   % y angle (Y)
phi2 = q(14);      % x angle (x) -->> ZYX pair yaw-pitch roll pair oder is imporant!!!
n12 = q(15); 
n22 = q(16);

% x_dot2 = 0.2;
% y_dot2 = 0.2;
% z_dot2 = 0;
% phi_dot2 = 0.4;   % \dotXangle
% theta_dot2 = 0;    % \dotYangle
% psi_dot2 = 0;      % \dotZangle in the order ZYX [psi2 theta2 phi2] -->> [Zangle Yangle Xangle] -->> Yaw , Pitch Roll 
% n1_dot1 = 0;
% n2_dot1 = 0;

mUAV = 5.0; 
mLink1 = 0.5; 
mLink2 = 0.5;
% l1 = 0.02; 
% l2 = 0.02;
%  
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
% p_b_dot2 = [x_dot2 y_dot2 z_dot2]';
% n_dot2 = [n1_dot1 n2_dot1]'; % joint angles rate 
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

% position of end effector 
p2_b2 = P102 + R012*p122 + R022*p2T2; % in base frame
% disp(p2_b1);
% disp(x1);
% disp(y1);
% disp(z1);
x2 = p_b2 + R_b2*p2_b2;  % in inertial frame
% disp('x2'); 
% disp(x2);


% J2c = [h1 h2; hat(h1)*(R01*p12) hat(h2)*p2c_b];
% Jacobian for end effector 
JT2 = [h12 h22; hat(h12)*R012*p122 hat(h22)*p2_b2];
JwT2 = JT2(1:3,1:2);
JvT2 = JT2(4:6,1:2);

% w_b2 = [phi_dot2*cos(psi2)*cos(theta2) - theta_dot2*sin(psi2); 
%       theta_dot2*cos(psi2) + phi_dot2*cos(theta2)*sin(psi2);
%                      psi_dot2 - phi_dot2*sin(theta2)];
% x2Dot = p_b_dot2 + hat(w_b2)*R_b2*p2_b12 + R_b2*JvT12*n_dot2;    % end effector

%% Load 
Mc = 0.5; 
e = [0 0 1]';
g = 9.81;
xc = [q(17);q(18);q(19)];
% xcDot = [q(10);q(11);q(12)];

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
% k = 100;
k = 2.5*10e4;
% f1 = k*norm(z1)*z1;
% f2 = k*norm(z2)*z2;

f1 = k*z1;
f2 = k*z2;

% f_sg = 0.5*Mc*g + 10;
McHat1= q(23);
McHat2 = q(24);
McHat = McHat1+McHat2;
f1d = [0 2.0 1.0*McHat1*g]';
f2d = [0 -2.0 1.0*McHat2*g]';

v1d = [0.2 0 0]';
v2d = [0.2 0 0]';

% gamma = 30*eye(3,3);
% F1 = -gamma*(x1Dot - v1d) + f1d;
% F2 = -gamma*(x2Dot - v2d) + f2d; 
% 
% x1dDot = (F1 - f1)/m1; 
% x2dDot = (F2 - f2)/m2;
% f_z = [0 0 20]';
% Mc = q(23);
xcDot = (f1 + f2 - Mc*e*g)/Mc;
% xcDot = (f1 + f2 - Mc*e*g)/Mc;
% 
K1 = [.1 0 0;
      0 .1 0;
      0 0 .1];
% K1 = 0.1;
K2 = K1;
x1Dot = K1*(f1d - f1) + v1d;
% x1Dot
x2Dot = K2*(f2d - f2) + v2d;
% x2Dot

% Km = [1 0 0;
%       0 1 0; 
%       0 0 10];
% % McHat1= q(23);
% McHat2 = q(24);
% McHat = McHat1+McHat2;
% McHatDot1 = -0.5*((x1Dot-v1d)'*g*e);
% McHatDot2 = -0.5*((x2Dot-v2d)'*g*e);

y_k =  1.0;
McHatDot1 = -0.6*((x1Dot-v1d)'*g*e);
McHatDot2 = -0.4*((x2Dot-v2d)'*g*e);

J_i1 = [eye(3,3),-hat(R_b1*p2_b1)*T_b1,R_b1*JvT1];
Ja1 = J_i1(:,1); 
Jb1 = J_i1(:,2);
Jf1 = J_i1(:,6);

J_New1 = [Ja1,Jb1,Jf1];

q1DotPrtl = inv(J_New1)*x1Dot;
q1Dot = [q1DotPrtl(1),q1DotPrtl(2),0,0,0,q1DotPrtl(3),0,0]';

J_i2 = [eye(3,3),-hat(R_b2*p2_b2)*T_b2,R_b2*JvT2];
Ja2 = J_i2(:,1); 
Jb2 = J_i2(:,2);
Jf2 = J_i2(:,6);   % phiDot to phi -->> X angle which is 6th variable in generalized coordinates

J_New2 = [Ja2,Jb2,Jf2];

q2DotPrtl = inv(J_New2)*x2Dot;
q2Dot = [q2DotPrtl(1),q2DotPrtl(2),0,0,0,q2DotPrtl(3),0,0]';

zDot(1:8,1)= q1Dot;
zDot(9:16,1) = q2Dot;
zDot(17:19,1)= q(20:22,1);
zDot(20:22,1)= xcDot;
zDot(23) = McHatDot1; 
zDot(24) = McHatDot2; 


