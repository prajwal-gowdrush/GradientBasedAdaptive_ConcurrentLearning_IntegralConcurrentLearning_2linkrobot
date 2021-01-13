function Xdot=ode_adaptintconlearn(t,X)
global p1 p2 p3 fd1 fd2 K1 K2 alpha gamma intofY4 Ts YsYs_sum YsTs_sum flagintcon timehistory tauhistory qhistory qdothistory del_t prevcheck  tvec mineigvec

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

%Adaptive Update Law
if flagintcon==1
 thetahatdot=gamma*Y2'*r + K2*gamma*(YsTs_sum - YsYs_sum*thetahat);  
else
 thetahatdot=gamma*Y2'*r;
end

%Matrices appearing in the dynamics
M= [ p1+2*p3*cos(q2), p2+p3*cos(q2); p2+p3*cos(q2),  p2]; %Inertia Matrix
V= [ -p3*sin(q2)*q2dot, -p3*sin(q2)*(q1dot+q2dot); p3*sin(q2)*q1dot , 0]; %Centrifugal Coriolis Matrix
Fd=[fd1 , 0;  0 , fd2];

%% Dynamics
Xdot(1,1)= q1dot;
Xdot(2,1)= q2dot;
Xdot(3:4,1)= M\(tau - V*qdot - Fd*qdot); % M\ is inverse of M
Xdot(5:9,1)=thetahatdot;

%Using the time history of states for the integral concurrent learning part of the adaptive update law
if (t-prevcheck)<= del_t
    if ~isempty(timehistory)
         if (t-timehistory(end))>10^-10
          intofY4 = intofY4  +  Y4(qhistory(end,:),qdothistory(end,:))*(t-timehistory(end));
          Ts = Ts + tauhistory(end,:)*(t-timehistory(end,:));
          
          timehistory=[timehistory;t];
          qhistory=[qhistory;q'];
          qdothistory=[qdothistory; qdot'];
          tauhistory=[tauhistory;tau'];
         end
    else
    timehistory=t;
    qhistory=q';
    qdothistory=qdot';
    tauhistory=tau';
    end   
else
Y3 = Y3_partial(qhistory(end,:),qdothistory(end,:)) - Y3_partial(qhistory(1,:),qdothistory(1,:));
Ys = Y3 + intofY4;
        if isempty(YsYs_sum)
           YsYs_sum= Ys'*Ys;
           YsTs_sum= Ys'*Ts';
           tvec=t;
           mineigvec=min(eig(YsYs_sum));
        elseif min(eig(YsYs_sum+ Ys'*Ys)) > min(eig(YsYs_sum))
           YsYs_sum= YsYs_sum + Ys'*Ys;
           YsTs_sum= YsTs_sum + Ys'*Ts';
           tvec=[tvec;t];
           mineigvec=[mineigvec;min(eig(YsYs_sum))];
           lambda=0.5;
            if min(eig(YsYs_sum))> lambda
               flagintcon=1;
            end
        end
prevcheck=t;
timehistory=[];
tauhistory=[];
qhistory=[];
qdothistory=[];
intofY4=0;
Ts=0;
end
