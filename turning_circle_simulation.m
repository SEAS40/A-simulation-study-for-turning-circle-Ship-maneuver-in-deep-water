
%Ship maneuver in deep water: A simulation study for turning circle
%Rohit Suryawanshi

clear; close all; clc;

%inital data
t_initial = 0;
u = 20/1.94384;   %conversion fron knots to m/s
v = 0;
r = 0;
psi = 0;
x = 0;
y = 0;

% calling the functions defined at bottom part of this code and plotting
for i = 1:1500
    [uI(i), vI(i), rI(i), psiI(i), xI(i), yI(i), timeI(i), accI(:,i),deltaI(i)] = equation_of_motion_deep(u,v,r,psi,x,y,t_initial);
    u = uI(i);
    v = vI(i);
    r = rI(i);
    psi = psiI(i);
    x = xI(i);
    y = yI(i);
    t_initial = timeI(i);
    acc = accI(:,i);
    d=deltaI(i);
end

%Finding particular points on turning circle to find turning circle
%particulars i.e. advance etc

tadv=find(psiI<-1.57 & psiI>-1.58);  % for 90 deg change in heading
xadv=xI(tadv)
yadv=yI(tadv)
xmax=max(xI)
ymax=max(yI)
tactical_rad=abs(xmax-0)
ttac=find(psiI>-3.15 & psiI<-3.14)   % for 180 deg change in heading
xtac=xI(ttac)
ytac=yI(ttac)

% loading refernce data for compariosn
data=xlsread('ferry_turning_circle.csv');  
realx=data(:,1);
realy=data(:,2);
figure(1)

% plotting the results as identified by the titles
plot(-yI,xI,'k','LineWidth',1.2)
grid on
hold on 
plot(-yI(tadv), xI(tadv), '*r')
hold on 
plot(-yI(ttac), xI(ttac), '*b')
xlabel('yo[m]'),ylabel('xo[m]')
title('Turning circle in deep water')
hold on;
scatter(realx,realy, 'r', 'LineWidth', 1.2,'SizeData', 0.8);
lgd=legend('Turning circle deep water','90 deg change of heading','180 deg change of heading','Turning circle validation data')
axis ([-100 1000 -300 800])
axis equal                      
hold off;
grid on

figure(2)
subplot(3,1,1);
plot(timeI,uI,'k','linewidth',1.2);
title('Longitudinal velocity')
xlabel('Time [s]'),ylabel('Velocity [m/s]')
legend ('u')
grid on
subplot(3,1,2);
plot(timeI,vI,'b','linewidth',1.2);
title('Transverse velocity')
xlabel('Time [s]'),ylabel('Velocity [m/s]')
legend ('v')
grid on
subplot(3,1,3);
plot(timeI,rI,'r','linewidth',1.2);
title('Turning rate')
xlabel('Time [s]'),ylabel('Turning rate [rad/s]')
grid on
legend ('r')
grid on
figure (3)
subplot(3,1,1);
plot(timeI,accI(1,:),'k','linewidth',1.2);
title('Longitudinal acceleration')
xlabel('Time [s]'),ylabel('Acceleration [m^2/s]')
grid on
legend ('up')
grid on
subplot(3,1,2);
plot(timeI,accI(2,:),'b','linewidth',1.2);
title('Transverse acceleration')
xlabel('Time (s)'),ylabel('Acceleration (m^2/s)')
grid on
legend ('vp')
grid on
subplot(3,1,3);
plot(timeI,accI(3,:),'r','linewidth',1.2);
title('Yaw acceleration')
xlabel('Time [s]'),ylabel('Acceleration [rad^2/s]')
grid on
legend ('rp')
grid on

figure(4)
plot(timeI,radtodeg(deltaI),'k','linewidth',1.2);
xlim([0 600]);
ylim([0 50]);
title('Rudder angle variation with time')
xlabel('Time [s]'),ylabel('Angle[deg]')
grid on
legend ('Rudder Angle')
grid on



function [u,v,r,psi,x,y,t,acceleration,d] = equation_of_motion_deep(u,v,r,psi,x,y,t)
rudder_rate = (2.3*pi)/180;
U0 = 20/1.94384;            %Initial velocity
L = 16*8.725;               % Ship length
U = sqrt(u^2 + v^2);        % combined velocity
t_final = 360;              % end time in seconds
N = 1000;                   % number of steps
dt = t_final/N;             

if rudder_rate*t <= ((35*pi)/180) %Defining rudder angle variation function
    d = rudder_rate*t;
else
    d = ((35*pi)/180);
end

%Finding longitundal and tranverse velocity in earth fixed coordinate
xp = u*cos(psi)-v*sin(psi);
yp = u*sin(psi)+v*cos(psi);

%Finding coordinate for ship position in earth fixed coordinate system
x = x+dt*xp;
y = y+dt*yp;
psi = psi+dt*r;

%Making dimension less
u = u/U;
v = v/U;
r = r*L/U;
du = u-U0/U;

%Some Ship data dimensionless
m  = 6765*10^-6;
I = 319*10^-6;
xg = (-116*10^-6)/m;


%Hydrodynamic coefficients & forces
Xu    = -4336*10^-6;
Xuu   = -2355*10^-6;
Xuuu  = -2594*10^-6;
Xvv   = -3279*10^-6;
Xr    = -19*10^-6;
Xrr   = -571*10^-6;
Xdd   = -2879*10^-6;
Xdddd = 2185*10^-6;
Xvvu  = -2559*10^-6;
Xrru  = -734*10^-6;
Xddu  = 3425*10^-6;
Xvr   = 4627*10^-6;
Xvd   = 877*10^-6;
Xrd   = -351*10^-6;
X = Xu.*du + Xuu.*du.^2 + Xuuu.*du.^3 + Xvv.*v.^2 + Xr.*r + Xrr.*r.^2 + Xdd.*d.^2 + Xdddd.*d.^4 + Xvvu.*du.*v.^2 + Xrru.*du.*r.^2 + Xddu.*du.*d.^2 + Xvr.*v.*r + Xvd.*v.*d + Xrd.*r.*d;
Yu    = 57*10^-6;
Yv    = -12095*10^-6;
Yvvv  = -137302*10^-6;
Yvp   = -7396*10^-6;
Yr    = 1901*10^-6;
Yrrr  = -1361*10^-6;
Yrp   = -600*10^-6;
Yd    = 3587*10^-6;
Ydd   = 98*10^-6;
Yddddd= -6262*10^-6;
Yru   = -1297*10^-6;
Ydu   = -5096*10^-6;
Ydddu = 3192*10^-6;
Yvrr  = -44365*10^-6;
Yvvr  = -36490*10^-6;
Yvdd  = 2199*10^-6;
Yrdd  = -2752*10^-6;
Y = Yu.*du + Yv.*v + Yvvv.*v.^3 + Yr.*r + Yrrr.*r.^3 + Yd.*d + Ydd.*d.^2 + Yddddd.*d.^5 + Yru.*r.*du + Ydu.*du.*d + Ydddu.*du.*d.^3 + Yvrr.*v.*r.^2 + Yvvr.*r.*v.^2 + Yvdd.*v.*d.^2 + Yrdd.*r.*d.^2;
Nu    = -36*10^-6;
Nv    = -3919*10^-6;
Nvvv  = -33857*10^-6;
Nvp   = 426*10^-6;
Nvpvv = 10049*10^-6;
Nr    = -2579*10^-6;
Nrrr  = -2253*10^-6;
Nrp   = -231*10^-6;
Nd    = -1621*10^-6;
Ndd   = -73*10^-6;
Nddddd= 2886*10^-6;
Nvu   = -3666*10^-6;
Nrrru = -1322*10^-6;
Ndu   = 2259*10^-6;
Ndddu = -1382*10^-6;
Nvvr  = -60110*10^-6;
Nvdd  = 570*10^-6;
Nvvd  = -2950*10^-6;
Nrdd  = 237*10^-6;
Nrrd  = -329*10^-6;

%Yawstabilitytest
lrdash = (Nr- m*xg)/(Yr-m);
lvdash=(Nv/Yv);
yaw=lrdash/lvdash

%Setting up the system of equation
N = Nu.*du + Nv.*v + Nvvv.*v.^3 + Nr.*r + Nrrr.*r.^3 + Nd.*d + Ndd.*d.^2 + Nddddd.*d.^5 + Nvu.*du.*v + Nrrru.*du.*r.^3 + Ndu.*d.*du + Ndddu.*du.*d.^3 + Nvvr.*r.*v.^2 + Nvdd.*v.*d.^2 + Nvvd.*d.*v.^2 + Nrdd.*r.*d.^2 + Nrrd.*d.*r.^2;
A = [m,0,0;0,m-Yvp,m*xg-Yrp;0,m*xg-Nvp-Nvpvv*v^2,I-Nrp];
b = [X+m*v*r+m*xg*r^2;Y-m*u*r;N-m*xg*u*r];
acceleration = A^-1*b;

%Finding velcoities using euler integration & Making results dimensional
u=u+dt*acceleration(1,1)*U/L;
v=v+dt*acceleration(2,1)*U/L;
r=r+dt*acceleration(3,1)*U/L;
u=u*U;
v=v*U;
r=r*U/L;
t=t+dt;
end

