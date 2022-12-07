% LMI for Controller Design of Multiple Solar PV Units Connected to Distribution Networks
Ri = [7.54*10^-3; 0.12; 0.12];
L = [34.24*10^-6;0.19*10^-3;1.05*10^-3];
Eg = 25*10^3;
fg = 60;
Vdci = 500;
w = 1;
Ci = 1;
A = [-sum(Ri)/sum(L) w 0; -w -sum(Ri)/sum(L) 0; 0 0 0];
B = [Vdci/sum(L) 0; 0 Vdci/sum(L); -1/Ci -1/Ci];
Q = sdpvar(3);
R = sdpvar(2);
P = sdpvar(3);
mat = [A'*P+P*A+Q, P*B; B'*P, R];
con = [P>=0,Q>=0,R>=0,mat<=0];
obj = [];
opt = sdpsettings('solver','sedumi');
sol = optimize(con,obj,opt);
K = -inv(value(R))*value(B)'*value(P)
