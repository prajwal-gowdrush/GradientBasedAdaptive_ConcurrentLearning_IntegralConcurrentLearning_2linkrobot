function Xdot=ode_stdadapt(t,X)
global K1 alpha gamma
%Assigning states from the input arguments to other variables
q1=X(1); q2=X(2); q1dot=X(3); q2dot=X(4); 
p1hat=X(5); p2hat=X(6); p3hat=X(7); fd1hat=X(8); fd2hat=X(9);
q=[q1;q2]; qdot=[q1dot;q2dot]; 
thetahat=[p1hat;p2hat;p3hat;fd1hat;fd2hat];

%Defining the desired trajectory
qd1=cos(0.5*t); qd2=2*cos(t);
qd1dot=-0.5*sin(0.5*t); qd2dot=-2*sin(t);
qd1doubledot=-0.25*cos(0.5*t); qd2doubledot=-2*cos(t);
qd=[qd1;qd2]; 
qd_dot=[qd1dot;qd2dot];
qd_doubledot=[qd1doubledot;qd2doubledot];

%Error definitions
e=q-qd; 
edot=qdot-qd_dot;
r= edot + alpha*e;

%% Regression Matrices: 
%Regression Matrix Y2
Y2=[alpha*(q1dot - qd1dot) - qd1doubledot, alpha*(q2dot - qd2dot) - qd2doubledot, 2*alpha*cos(q2)*(q1dot - qd1dot) - qd2doubledot*cos(q2) - 2*qd1doubledot*cos(q2) + alpha*cos(q2)*(q2dot - qd2dot) + sin(q2)*(qd2dot - alpha*(q2 - qd2))*(q1dot + q2dot) + q2dot*sin(q2)*(qd1dot - alpha*(q1 - qd1)), -q1dot,  0;    0, alpha*(q1dot - qd1dot) - qd2doubledot - qd1doubledot + alpha*(q2dot - qd2dot), alpha*cos(q2)*(q1dot - qd1dot) - qd1doubledot*cos(q2) - q1dot*sin(q2)*(qd1dot - alpha*(q1 - qd1)), 0, -q2dot];
 %Control Law
tau= -Y2*thetahat - e - K1*r; 

%Actual Parameter values (Have not been used in the control law)
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;
fd1 = 5.3;
fd2 = 1.1;


%Adaptive Update Law
thetahatdot=gamma*Y2'*r;

%Matrices appearing in the dynamics
M= [ p1+2*p3*cos(q2), p2+p3*cos(q2); p2+p3*cos(q2),  p2]; %Inertia Matrix
V= [ -p3*sin(q2)*q2dot, -p3*sin(q2)*(q1dot+q2dot); p3*sin(q2)*q1dot , 0]; %Centrifugal Coriolis Matrix
Fd=[fd1 , 0;  0 , fd2];

%% Dynamics
Xdot(1,1)= q1dot;
Xdot(2,1)= q2dot;
Xdot(3:4,1)= M\(tau - V*qdot - Fd*qdot); % M\ is inverse of M
Xdot(5:9,1)=thetahatdot;