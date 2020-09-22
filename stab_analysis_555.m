format long

s = dlmread('555-pi-1.dat');
a=s(:,1);
E=s(:,2);

syms p1 p2 p3 p4 p5 q0 q1 q2 q3 q4 q5 r0 r1 r2 r3 r4 r5

eqn1 = (E(1)^2)*a(1)*p1 + (E(1)^2*a(1)^2)*p2 + (E(1)^2*a(1)^3)*p3 + (E(1)^2*a(1)^4)*p4 + (E(1)^2*a(1)^5)*p5  + E(1)*q0 + E(1)*a(1)*q1 + E(1)*(a(1)^2)*q2 + E(1)*(a(1)^3)*q3 + E(1)*(a(1)^4)*q4 + E(1)*(a(1)^5)*q5  + r0 + a(1)*r1 + (a(1)^2)*r2 + (a(1)^3)*r3 + (a(1)^4)*r4 + (a(1)^5)*r5  == -(E(1)^2);

eqn2 = (E(2)^2)*a(2)*p1 + (E(2)^2*a(2)^2)*p2 + (E(2)^2*a(2)^3)*p3 + (E(2)^2*a(2)^4)*p4 + (E(2)^2*a(2)^5)*p5  + E(2)*q0 + E(2)*a(2)*q1 + E(2)*(a(2)^2)*q2 + E(2)*(a(2)^3)*q3 + E(2)*(a(2)^4)*q4 + E(2)*(a(2)^5)*q5  + r0 + a(2)*r1 + (a(2)^2)*r2 + (a(2)^3)*r3 + (a(2)^4)*r4 + (a(2)^5)*r5  == -(E(2)^2);

eqn3 = (E(3)^2)*a(3)*p1 + (E(3)^2*a(3)^2)*p2 + (E(3)^2*a(3)^3)*p3 + (E(3)^2*a(3)^4)*p4 + (E(3)^2*a(3)^5)*p5  + E(3)*q0 + E(3)*a(3)*q1 + E(3)*(a(3)^2)*q2 + E(3)*(a(3)^3)*q3 + E(3)*(a(3)^4)*q4 + E(3)*(a(3)^5)*q5  + r0 + a(3)*r1 + (a(3)^2)*r2 + (a(3)^3)*r3 + (a(3)^4)*r4 + (a(3)^5)*r5  == -(E(3)^2);

eqn4 = (E(4)^2)*a(4)*p1 + (E(4)^2*a(4)^2)*p2 + (E(4)^2*a(4)^3)*p3 + (E(4)^2*a(4)^4)*p4 + (E(4)^2*a(4)^5)*p5  + E(4)*q0 + E(4)*a(4)*q1 + E(4)*(a(4)^2)*q2 + E(4)*(a(4)^3)*q3 + E(4)*(a(4)^4)*q4 + E(4)*(a(4)^5)*q5  + r0 + a(4)*r1 + (a(4)^2)*r2 + (a(4)^3)*r3 + (a(4)^4)*r4 + (a(4)^5)*r5  == -(E(4)^2);

eqn5 = (E(5)^2)*a(5)*p1 + (E(5)^2*a(5)^2)*p2 + (E(5)^2*a(5)^3)*p3 + (E(5)^2*a(5)^4)*p4 + (E(5)^2*a(5)^5)*p5  + E(5)*q0 + E(5)*a(5)*q1 + E(5)*(a(5)^2)*q2 + E(5)*(a(5)^3)*q3 + E(5)*(a(5)^4)*q4 + E(5)*(a(5)^5)*q5  + r0 + a(5)*r1 + (a(5)^2)*r2 + (a(5)^3)*r3 + (a(5)^4)*r4 + (a(5)^5)*r5  == -(E(5)^2);

eqn6 = (E(6)^2)*a(6)*p1 + (E(6)^2*a(6)^2)*p2 + (E(6)^2*a(6)^3)*p3 + (E(6)^2*a(6)^4)*p4 + (E(6)^2*a(6)^5)*p5  + E(6)*q0 + E(6)*a(6)*q1 + E(6)*(a(6)^2)*q2 + E(6)*(a(6)^3)*q3 + E(6)*(a(6)^4)*q4 + E(6)*(a(6)^5)*q5  + r0 + a(6)*r1 + (a(6)^2)*r2 + (a(6)^3)*r3 + (a(6)^4)*r4 + (a(6)^5)*r5  == -(E(6)^2);

eqn7 = (E(7)^2)*a(7)*p1 + (E(7)^2*a(7)^2)*p2 + (E(7)^2*a(7)^3)*p3 + (E(7)^2*a(7)^4)*p4 + (E(7)^2*a(7)^5)*p5  + E(7)*q0 + E(7)*a(7)*q1 + E(7)*(a(7)^2)*q2 + E(7)*(a(7)^3)*q3 + E(7)*(a(7)^4)*q4 + E(7)*(a(7)^5)*q5  + r0 + a(7)*r1 + (a(7)^2)*r2 + (a(7)^3)*r3 + (a(7)^4)*r4 + (a(7)^5)*r5  == -(E(7)^2);

eqn8 = (E(8)^2)*a(8)*p1 + (E(8)^2*a(8)^2)*p2 + (E(8)^2*a(8)^3)*p3 + (E(8)^2*a(8)^4)*p4 + (E(8)^2*a(8)^5)*p5  + E(8)*q0 + E(8)*a(8)*q1 + E(8)*(a(8)^2)*q2 + E(8)*(a(8)^3)*q3 + E(8)*(a(8)^4)*q4 + E(8)*(a(8)^5)*q5  + r0 + a(8)*r1 + (a(8)^2)*r2 + (a(8)^3)*r3 + (a(8)^4)*r4 + (a(8)^5)*r5  == -(E(8)^2);

eqn9 = (E(9)^2)*a(9)*p1 + (E(9)^2*a(9)^2)*p2 + (E(9)^2*a(9)^3)*p3 + (E(9)^2*a(9)^4)*p4 + (E(9)^2*a(9)^5)*p5  + E(9)*q0 + E(9)*a(9)*q1 + E(9)*(a(9)^2)*q2 + E(9)*(a(9)^3)*q3 + E(9)*(a(9)^4)*q4 + E(9)*(a(9)^5)*q5  + r0 + a(9)*r1 + (a(9)^2)*r2 + (a(9)^3)*r3 + (a(9)^4)*r4 + (a(9)^5)*r5  == -(E(9)^2);

eqn10 = (E(10)^2)*a(10)*p1 + (E(10)^2*a(10)^2)*p2 + (E(10)^2*a(10)^3)*p3 + (E(10)^2*a(10)^4)*p4 + (E(10)^2*a(10)^5)*p5  + E(10)*q0 + E(10)*a(10)*q1 + E(10)*(a(10)^2)*q2 + E(10)*(a(10)^3)*q3 + E(10)*(a(10)^4)*q4 + E(10)*(a(10)^5)*q5  + r0 + a(10)*r1 + (a(10)^2)*r2 + (a(10)^3)*r3 + (a(10)^4)*r4 + (a(10)^5)*r5  == -(E(10)^2);

eqn11 = (E(11)^2)*a(11)*p1 + (E(11)^2*a(11)^2)*p2 + (E(11)^2*a(11)^3)*p3 + (E(11)^2*a(11)^4)*p4 + (E(11)^2*a(11)^5)*p5  + E(11)*q0 + E(11)*a(11)*q1 + E(11)*(a(11)^2)*q2 + E(11)*(a(11)^3)*q3 + E(11)*(a(11)^4)*q4 + E(11)*(a(11)^5)*q5  + r0 + a(11)*r1 + (a(11)^2)*r2 + (a(11)^3)*r3 + (a(11)^4)*r4 + (a(11)^5)*r5  == -(E(11)^2);

eqn12 = (E(12)^2)*a(12)*p1 + (E(12)^2*a(12)^2)*p2 + (E(12)^2*a(12)^3)*p3 + (E(12)^2*a(12)^4)*p4 + (E(12)^2*a(12)^5)*p5  + E(12)*q0 + E(12)*a(12)*q1 + E(12)*(a(12)^2)*q2 + E(12)*(a(12)^3)*q3 + E(12)*(a(12)^4)*q4 + E(12)*(a(12)^5)*q5  + r0 + a(12)*r1 + (a(12)^2)*r2 + (a(12)^3)*r3 + (a(12)^4)*r4 + (a(12)^5)*r5  == -(E(12)^2);

eqn13 = (E(13)^2)*a(13)*p1 + (E(13)^2*a(13)^2)*p2 + (E(13)^2*a(13)^3)*p3 + (E(13)^2*a(13)^4)*p4 + (E(13)^2*a(13)^5)*p5  + E(13)*q0 + E(13)*a(13)*q1 + E(13)*(a(13)^2)*q2 + E(13)*(a(13)^3)*q3 + E(13)*(a(13)^4)*q4 + E(13)*(a(13)^5)*q5  + r0 + a(13)*r1 + (a(13)^2)*r2 + (a(13)^3)*r3 + (a(13)^4)*r4 + (a(13)^5)*r5  == -(E(13)^2);

eqn14 = (E(14)^2)*a(14)*p1 + (E(14)^2*a(14)^2)*p2 + (E(14)^2*a(14)^3)*p3 + (E(14)^2*a(14)^4)*p4 + (E(14)^2*a(14)^5)*p5  + E(14)*q0 + E(14)*a(14)*q1 + E(14)*(a(14)^2)*q2 + E(14)*(a(14)^3)*q3 + E(14)*(a(14)^4)*q4 + E(14)*(a(14)^5)*q5  + r0 + a(14)*r1 + (a(14)^2)*r2 + (a(14)^3)*r3 + (a(14)^4)*r4 + (a(14)^5)*r5  == -(E(14)^2);

eqn15 = (E(15)^2)*a(15)*p1 + (E(15)^2*a(15)^2)*p2 + (E(15)^2*a(15)^3)*p3 + (E(15)^2*a(15)^4)*p4 + (E(15)^2*a(15)^5)*p5  + E(15)*q0 + E(15)*a(15)*q1 + E(15)*(a(15)^2)*q2 + E(15)*(a(15)^3)*q3 + E(15)*(a(15)^4)*q4 + E(15)*(a(15)^5)*q5  + r0 + a(15)*r1 + (a(15)^2)*r2 + (a(15)^3)*r3 + (a(15)^4)*r4 + (a(15)^5)*r5  == -(E(15)^2);

eqn16 = (E(16)^2)*a(16)*p1 + (E(16)^2*a(16)^2)*p2 + (E(16)^2*a(16)^3)*p3 + (E(16)^2*a(16)^4)*p4 + (E(16)^2*a(16)^5)*p5  + E(16)*q0 + E(16)*a(16)*q1 + E(16)*(a(16)^2)*q2 + E(16)*(a(16)^3)*q3 + E(16)*(a(16)^4)*q4 + E(16)*(a(16)^5)*q5  + r0 + a(16)*r1 + (a(16)^2)*r2 + (a(16)^3)*r3 + (a(16)^4)*r4 + (a(16)^5)*r5  == -(E(16)^2);

eqn17 = (E(17)^2)*a(17)*p1 + (E(17)^2*a(17)^2)*p2 + (E(17)^2*a(17)^3)*p3 + (E(17)^2*a(17)^4)*p4 + (E(17)^2*a(17)^5)*p5  + E(17)*q0 + E(17)*a(17)*q1 + E(17)*(a(17)^2)*q2 + E(17)*(a(17)^3)*q3 + E(17)*(a(17)^4)*q4 + E(17)*(a(17)^5)*q5  + r0 + a(17)*r1 + (a(17)^2)*r2 + (a(17)^3)*r3 + (a(17)^4)*r4 + (a(17)^5)*r5  == -(E(17)^2);

[M,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10, eqn11, eqn12, eqn13, eqn14, eqn15, eqn16, eqn17], [p1, p2, p3, p4, p5, q0, q1, q2, q3, q4, q5, r0, r1, r2, r3, r4, r5]);

X = vpa(linsolve(M,B),8);

p1=X(1);
p2=X(2);
p3=X(3);
p4=X(4);
p5=X(5);
q0=X(6);
q1=X(7);
q2=X(8);
q3=X(9);
q4=X(10);
q5=X(11);
r0=X(12);
r1=X(13);
r2=X(14);
r3=X(15);
r4=X(16);
r5=X(17);

syms a

P = 1 + (a*p1) + (a^2*p2) + (a^3*p3) + (a^4*p4) + (a^5*p5);
Q = q0 + (a*q1) + (a^2*q2) + (a^3*q3) + (a^4*q4) + (a^5*q5);
R = r0 + (a*r1) + (a^2*r2) + (a^3*r3) + (a^4*r4) + (a^5*r5);

bp = (Q^2)-(4*P*R);

f = -(Q/(2*P))+(sqrt((Q^2)-(4*P*R)))/(2*P);

d=diff(f);
