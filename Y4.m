function y4=Y4(q,qdot)

q2=q(2); q1dot=qdot(1); q2dot=qdot(2);
y4=[ 0, 0, q2dot^2*sin(q2) + q1dot*q2dot*sin(q2) - q2dot*sin(q2)*(q1dot + q2dot), q1dot, 0; 0, 0, sin(q2)*q1dot^2 + q2dot*sin(q2)*q1dot, 0, q2dot];