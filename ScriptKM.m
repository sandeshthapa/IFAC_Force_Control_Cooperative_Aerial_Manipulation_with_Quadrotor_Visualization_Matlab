close all; clear all; clc; 
sim('KinematicModelG')

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
plot(t,f1(:,1),'--b',t,f2(:,1),':r',t,f1(:,2),'r',t,f2(:,2),'b',t,f1(:,3),'--k',t,f2(:,3),':b','LineWidth',1.5)
xlabel('time (sec)')
ylabel('Squeeze force (N)')
legend('f1x','f2x','f1y','f2y','f1z','f2z','Location','best','Orientation','horizontal');
hold on 

%%


figure
plot(t,q1(:,1),'r',t,q2(:,1),'b','LineWidth',2)
% plot(t,q(1))
xlabel('time (sec)')
ylabel('x postion (m)')
legend('Agent 1(LAM)','Agent 2(RAM)','Location','best');
hold on 

figure
plot(t,q1(:,2),'r',t,q2(:,2),'b','LineWidth',2)
xlabel('time (sec)')
ylabel('y postion (m)')
legend('Agent 1(LAM)','Agent 2(RAM)','Location','best');
hold on 

figure
plot(t,q1(:,6),'r',t,q2(:,6),'b','LineWidth',2)
xlabel('time (sec)')
ylabel('roll angle(ZYX X anlge) (rad))')
legend('Agent 1(LAM)','Agent 2(RAM)','Location','best');
hold on 

figure
plot(t,xc(:,1),'r',t,xc(:,2),'m',t,xc(:,3),'b','LineWidth',2)
xlabel('time (sec)')
ylabel('position of load (m)')
legend('x pos','y pos','z pos','Location','best');
hold on 

figure
plot(t,xcDot(:,1),'r',t,xcDot(:,2),'m',t,xcDot(:,3),'b','LineWidth',2)
% plot(t,q(:,10))
xlabel('time (sec)')
ylabel('velocity of load (m/s)')
legend('X velocity','Y Velocity','Z Velocity','Location','best');
hold on 

figure
plot(t,f1(:,1),'r',t,f2(:,1),'b','LineWidth',2)
xlabel('time (sec)')
ylabel('x direction squeeze force (N)')
legend('f1x','f2x','Location','best');
hold on 

figure
plot(t,f1(:,2),'r',t,f2(:,2),'b','LineWidth',2)
xlabel('time (sec)')
ylabel('y direction squeeze force (N)')
legend('f1y','f2y','Location','best');
hold on

figure
plot(t,f1(:,3),'r',t,f2(:,3),'b','LineWidth',2)
xlabel('time (sec)')
ylabel('z direction squeeze force (N)')
legend('f1z','f2z','Location','best');
hold on 

figure
plot(t,f1(:,1),'c',t,f2(:,1),'g',t,f1(:,2),'r',t,f2(:,2),'b',t,f1(:,3),'y',t,f2(:,3),'k','LineWidth',2)
xlabel('time (sec)')
ylabel('squeeze force (N)')
legend('f1x','f2x','f1y','f2y','f1z','f2z','Location','best','Orientation','horizontal');
hold on 

figure
plot(t,M1Hat,'k', t, M2Hat,'b','LineWidth',2)
xlabel('time (sec)')
ylabel('Mass of Unknown Payload')
hold on 






