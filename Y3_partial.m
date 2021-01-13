function y3_partial=Y3_partial(q,qdot)

q2=q(2); q1dot=qdot(1); q2dot=qdot(2);
y3_partial=[ q1dot, q2dot, 2*q1dot*cos(q2) + q2dot*cos(q2), 0, 0; 0, q1dot + q2dot, q1dot*cos(q2), 0, 0];
