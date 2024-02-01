% UCLA MAE294 HW2
% written by Yuan Hung Lo

clear all
close all
%% process of solution
% get Ktw and Mtw
% K=N*delta*S*N^-1 ***for the wire and blades
% M=N*delta*In*N^-1 ***for the rigid body
% get eigen vectors and values from Mtw^-1 * Ktw
% w=sqrt(eigen values)


syms k
z=[0 0 0]';
E=68e9;
G=25e9;
p=2698.9;
%% Calculating N1 S1 for the blade
% Calculate N
L1=[0.09875,0.2,0.025]';
n11=[0 0 -1]';
n21=[-1 0 0]';
n31=[0 1 0]';

N1=[n11           n21           n31           z   z   z;
    cross(L1,n11) cross(L1,n21) cross(L1,n31) n11 n21 n31];
%N1

% calculate S
b=0.05;
h=0.0025;
l=0.15;

ix=b*h^3/12;
iy=h*b^3/12;
A=b*h;
J=(b*h^3)/3 * (1-192/pi^5 * (h/b) * symsum(1/(2*k-1)^5 * tanh((2*k-1)*pi*b/(2*h)), k, 1, 5)); 
% n = 1 3 5 7 9 to simplify calculation

S1=[l/(E*ix)     0             0            0             -l^2/(2*E*ix) 0;
    0            l/(E*iy)      0            l^2/(2*E*iy)   0            0;
    0            0             l/(G*J)      0              0            0;
    0            l^2/(2*E*iy)  0            l^3/(3*E*iy)   0            0;
   -l^2/(2*E*ix) 0             0            0              l^3/(3*E*ix) 0;
    0            0             0            0              0            l/(E*A)];


%S1
%% Calculate N2 S2 for the wire
L2=[0.00125,0.2,0.04875]';
n12=[0 0 -1]';
n22=[-1 0 0]';
n32=[0 1 0]';

N2=[n12           n22           n32           z   z   z;
    cross(L2,n12) cross(L2,n22) cross(L2,n32) n12 n22 n32];
%N2

b=0.0025;
l=0.15;

J=(b^4)/3 * (1-192/pi^5 * symsum(1/(2*k-1)^5 * tanh((2*k-1)*pi/2), k, 1, 5)); 
% n = 1 3 5 7 9 to simplify calculation

S2=[12*l/(E*b^4) 0             0            0            -6*l^2/(E*b^4) 0;
    0            12*l/(E*b^4)  0            6*l^2/(E*b^4) 0            0;
    0            0             l/(G*J)      0              0            0;
    0            6*l^2/(E*b^4) 0            4*l^3/(E*b^4)  0            0;
    -6*l^2/(E*b^4) 0           0            0              4*l^3/(E*b^4) 0;
    0            0             0            0              0            l/(E*b^2)];

%% sum for Ktw
delt=[zeros(3,3) eye(3);
    eye(3) zeros(3,3)];
K1=N1* delt /S1/N1;
K2=N2* delt /S2/N2;

Ktw=vpa(K1+K2);



%% Calculate M
V=0.05*0.1*0.05;
b=0.05;
l=0.1;
h=0.05;

n1m=[1 0 0]';
n2m=[0 1 0]';
n3m=[0 0 1]';
Lm=[0.05 0.225 0.025]';
Nm=[n1m           n2m           n3m           z   z   z;
    cross(Lm,n1m) cross(Lm,n2m) cross(Lm,n3m) n1m n2m n3m];

Ix=p*V/12*(b^2+h^2);
Iy=p*V/12*(l^2+h^2);
Iz=p*V/12*(b^2+l^2);
In=diag([Ix Iy Iz p*V p*V p*V]);

Mtw=Nm*delt*In/(Nm);


%% get eigens for w and modes
[EigenVectors, EigenValues]=eig(Mtw\Ktw);



w=sqrt(EigenValues);
%w
%EigenVectors

%% using EulerStiffnessMatrix function to cross verify results
% constraints=[0.09875,0.2,0.025,0,1,0,1,0,0,0.15,0.05,0.0025,E,G,1,0; ...
%             0.00125,0.2,0.04875,0,1,0,1,0,0,0.15,0.0025,0.0025,E,G,1,0;
%             1,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0];
% KK=EulerStiffnessMatrix(constraints);
% 
% [EigenVectorsc, EigenValuesc]=eig(Mtw\vpa(KK));
% 
% 
% wc=sqrt(EigenValuesc);
%wc
%EigenVectorsc
% pretty close results!

%% Display Final Solution
modeShapes=["Mode Shape 1";"Mode Shape 2";"Mode Shape 3";"Mode Shape 4";"Mode Shape 5";"Mode Shape 6"];

w1=vpa([EigenVectors(1,6);EigenVectors(1,5);EigenVectors(1,4);EigenVectors(1,3);EigenVectors(1,2);EigenVectors(1,1)],6);
w2=vpa([EigenVectors(2,6);EigenVectors(2,5);EigenVectors(2,4);EigenVectors(2,3);EigenVectors(2,2);EigenVectors(2,1)],6);
w3=vpa([EigenVectors(3,6);EigenVectors(3,5);EigenVectors(3,4);EigenVectors(3,3);EigenVectors(3,2);EigenVectors(3,1)],6);
v1=vpa([EigenVectors(4,6);EigenVectors(4,5);EigenVectors(4,4);EigenVectors(4,3);EigenVectors(4,2);EigenVectors(4,1)],6);
v2=vpa([EigenVectors(5,6);EigenVectors(5,5);EigenVectors(5,4);EigenVectors(5,3);EigenVectors(5,2);EigenVectors(5,1)],6);
v3=vpa([EigenVectors(6,6);EigenVectors(6,5);EigenVectors(6,4);EigenVectors(6,3);EigenVectors(6,2);EigenVectors(6,1)],6);
Frequencies=vpa([w(6,6);w(5,5);w(4,4);w(3,3);w(2,2);w(1,1)],6);

Answer=table(modeShapes,w1,w2,w3,v1,v2,v3,Frequencies)
% will implement for loops and sorting function in the future 
% to make this process more efficient
