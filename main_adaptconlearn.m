%% MAINSCRIPT STANDARD ADAPTIVE CONTROLLER WITH CONCURRENT LEARNING
clc; clear; close all;

global K1 K2 alpha gamma tprevious qdotprevious Y1Y1_sum Y1T1_sum flagcon prevcheck tvec mineigvec

prevcheck=-0.3;
flagcon=0;
tprevious=[];
qdotprevious=[];
Y1Y1_sum=[];
Y1T1_sum=[];

%%%%%%%%%%%%%%%%%%%Tunable Controller Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%
K1=9; %Constant multiplying 'r' in the control law
K2=0.01;  %Constant multiplying the concurrent learning part of the adaptive update law 
alpha=1; %Constant in the definition of the error signal 'r'
gamma=diag([300; 2; 20; 20 ;20]);%Learning rate matrix in the adaptive update law

%% Augmented Initial State:
% Initial State:
q10=3; q20=4.5; q1dot0=0; q2dot0=0;
x0=[q10;q20;q1dot0;q2dot0];

% Initial Values of the adaptive estimates: 
%(We pick these based on our expectations of the actual parameter values)
p1hat0=8; p2hat0=0.7; p3hat0=0.5; fd1hat0=8; fd2hat0=2; 
thetahat0=[p1hat0;p2hat0;p3hat0;fd1hat0;fd2hat0];

%Augmented Initial State
X0=[x0;thetahat0];

%% Simulation time
t_sim=40;

%% Obtaining the solution to the ode based on the control law specified in ode_compadaptgrad
tspan=[0 t_sim];
options = odeset('AbsTol',1e-3,'RelTol',1e-3);
[time,X_Sol]=ode113(@ode_adaptconlearn,tspan,X0,options);

%% Extracting the trajectories of the states and the estimates from the ODE solution matrix X_Sol:
q1=X_Sol(:,1);
q2=X_Sol(:,2);
q1dot=X_Sol(:,3);
q2dot=X_Sol(:,4);
p1hat=X_Sol(:,5);
p2hat=X_Sol(:,6);
p3hat=X_Sol(:,7);
fd1hat=X_Sol(:,8);
fd2hat=X_Sol(:,9);

%% Computing desired trajectories, errors and controls from the ode solution for the plots:

%Desired trajectories:
qd1=cos(0.5*time); qd2=2*cos(time);
qd1dot=-0.5*sin(0.5*time); qd2dot=-2*sin(time);
qd1doubledot=-0.25*cos(0.5*time); qd2doubledot=-2*cos(time);

%Actual Parameter values (Have not been used in the control law)
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;
fd1 = 5.3;
fd2 = 1.1;

%% Errors
e1=q1-qd1; e2=q2-qd2;
e1dot=q1dot-qd1dot; e2dot=q2dot-qd2dot;
p1error= p1-p1hat; p2error= p2-p2hat; p3error= p3-p3hat; fd1error= fd1-fd1hat; fd2error= fd2-fd2hat;
r1=e1dot+alpha*e1;
r2=e1dot+alpha*e2;

%% Controls
y2thetahat_1=p3hat.*(2.*alpha.*cos(q2).*(q1dot - qd1dot) - qd2doubledot.*cos(q2) - 2.*qd1doubledot.*cos(q2) + alpha.*cos(q2).*(q2dot - qd2dot) + sin(q2).*(qd2dot - alpha.*(q2 - qd2)).*(q1dot + q2dot) + q2dot.*sin(q2).*(qd1dot - alpha.*(q1 - qd1))) - p2hat.*(qd2doubledot - alpha.*(q2dot - qd2dot)) - fd1hat.*q1dot - p1hat.*(qd1doubledot - alpha.*(q1dot - qd1dot));
y2thetahat_2= - fd2hat.*q2dot - p2hat.*(qd1doubledot + qd2doubledot - alpha.*(q1dot - qd1dot) - alpha.*(q2dot - qd2dot)) - p3hat.*(qd1doubledot.*cos(q2) - alpha.*cos(q2).*(q1dot - qd1dot) + q1dot.*sin(q2).*(qd1dot - alpha.*(q1 - qd1)));
tau1= -y2thetahat_1 - e1 - K1*r1; 
tau2= -y2thetahat_2 - e2 - K1*r2; 
%% Plots
% figure;
subplot(2,2,1);
plot(time,q1,'r',time,q2,'b',time,qd1,'r--',time,qd2,'b--','LineWidth',1.5);
xlabel('Time');
ylabel('States');
title('Desired Vs Actual Trajectories');
xlim([0 time(end)]);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
leg1 = legend('$q_1$','$q_2$','$q_{d1}$','$q_{d2}$');
set(leg1,'Interpreter','latex');
hold on;

% figure;
subplot(2,2,2);
plot(time,e1,'r',time,e2,'b','LineWidth',1.5);
xlabel('Time');
ylabel('Tracking Errors');
title('Tracking Errors Vs Time');
xlim([0 time(end)]);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
leg3 = legend('$e_1$','$e_2$');
set(leg3,'Interpreter','latex');
hold on;

% figure;
subplot(2,2,3);
plot(time,p1hat,'r',time,p2hat,'b',time,p3hat,'k',time,fd1hat,'m',time,fd2hat,'g','LineWidth',1.5);
xlabel('Time');
ylabel('Adaptive Estimates');
title('Adaptive Estimates Vs Time');
xlim([0 time(end)]);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
leg2 = legend('$\widehat{p}_1$','$\widehat{p}_2$','$\widehat{p}_3$','$\widehat{f}_{d1}$','$\widehat{f}_{d2}$');
set(leg2,'Interpreter','latex');
hold on;

% figure;
subplot(2,2,4);
plot(time,p1error,'r',time,p2error,'b',time,p3error,'k',time,fd1error,'m',time,fd2error,'g','LineWidth',1.5);
xlabel('Time');
ylabel('Actual value - Adaptive Estimate');
title('Parameter Estimate Errors');
xlim([0 time(end)]);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
leg2 = legend('$p_1-\widehat{p}_1$','$p_2-\widehat{p}_2$','$p_3-\widehat{p}_3$','$f_{d1}-\widehat{f}_{d1}$','$f_{d2}-\widehat{f}_{d2}$');
set(leg2,'Interpreter','latex');
hold on;


figure;
plot(time,tau1,'r',time,tau2,'b','LineWidth',1.5);
xlabel('Time');
ylabel('Controls');
title('Controls Vs Time');
xlim([0 time(end)]);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
leg4 = legend('$\tau_1$','$\tau_2$');
set(leg4,'Interpreter','latex');
hold on;

figure;
plot(tvec,mineigvec,'b','LineWidth',1.5);
xlabel('Time');
ylabel('Minimum Eigenvalue of Y''Y summation');
title('\lambda_{min} Vs Time');
xlim([0 time(end)]);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
hold on;

% The following lines are useful for tuning the gains:
maxlink1torque=max(abs(tau1))
maxlink2torque=max(abs(tau2))