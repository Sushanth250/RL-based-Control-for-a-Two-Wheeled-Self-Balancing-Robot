%% Multi Step Q Learning for Quanser Helicopter
clc;
clear all;
close all;

m = 1.3872;
l_cm = 0.1860;
Jp= 0.01778;
Jy = 0.0084;
Bp=0.8;
By = 0.3180;
Kpp = 0.2040;
Kyy = 0.0720;
Kpy = 0.0068;
Kyp = 0.0219;

Jtp=Jp + (m*l_cm*l_cm);
Jty = Jy + (m*l_cm*l_cm);

% Define the system matrices
A = [0 0 1 0;
     0 0 0 1;
     0 0 -Bp/(Jtp) 0;
     0 0 0 -By/(Jty)]

B = [0 0;
     0 0;
     Kpp/Jtp Kpy/Jtp;
     Kyp/Jty Kyy/Jty];

C = [1 0 0 0;
     0 1 0 0];
D = 0;

% Create state-space model object
sys = ss(A, B, C, D);

% Discretize the system
Ts = 0.1;
sysd = c2d(sys, Ts,'zoh');

% Initialize Q-table
S = [1 0;
     0 1];

R = [1 0;
     0 1];
sysd.C
sysd.A
sysd.B
Q = sysd.C' * S * sysd.C
%% Finding Solution based on ARE for the H Matrix

% Solution of the ARE (Algebric Riccati Equation)
P=idare(sysd.A,sysd.B,Q,R);
P;

% Solution Based on ARE for the H (G in Paper of MsQL) Matrix
G_optimal=[Q+sysd.A'*P*sysd.A sysd.A'*P*sysd.B; (sysd.A'*P*sysd.B)' R+sysd.B'*P*sysd.B]	
k=-inv(R+sysd.B'*P*sysd.B)*(sysd.A'*P*sysd.B)'

% K = [1 1 1 1;
%      1 1 1 1];

% G11 generated randomly form numbers btw [0,1]
G11 = rand(4,4) ;

G12 = rand(4,2);
G22 = rand(2,2);

K =  -inv(G22)*G12';

G_d = [G11 G12;
       G12' G22];
G_d


eta = rand(1,50);
rho = rand(28,50);
x(1:4,1:400) = 0;
for p = 1 : 300
    x(:,p) = [25 -40 0 0]';
end
%x(:,1) = [25 -40 0 0]';
 

%% Multi Step Q Learning

% Size of the Dataset
M = 50;

% Step Size (Ni)
% Ni = 1+4*sqrt(i) for the ith iteration


% Run the loop for 300 iterations
for i = 1:300
    a1 = 1;
    noise1= a1*0.01*(0.1*sin(1*i)^2*cos(9*i)+sin(2*i)^2*cos(0.1*i)+sin(-1.2*i)^2*cos(0.5*i)+sin(i)^5+sin(1.12*i)^2+cos(2.4*i)*sin(2.4*i)^3);
    N_i =floor( 1 + (4*sqrt(i)) );
    % Step - 2
    u=K*x(:,i) + noise1;


    % The next State
    x(:,i+1)=sysd.A*x(:,i)+sysd.B*u;

    % Step - 3 (Collecting Measured Dataset- S50)
    
    for k=1:49
        eta(1,k) = eta(1,k+1);
    end

    for k=1:49
        rho(:,k) = rho(:,k+1);
    end
    
    %rho 50
    kron( x(:,i)',x(:,i)');
    size(kron( x(:,i)',x(:,i)'));
    size((K*x(:,i))');
    size(kron((K*x(:,i))' , (K*x(:,i))'));
     size(2*kron( x(:,i)',(K*x(:,i))'));
     
    temp2 = [kron( x(:,i)',x(:,i)'),  2*kron( x(:,i)',(K*x(:,i))'), kron((K*x(:,i))' , (K*x(:,i))')];
    rho(:,50) = temp2;
    
    size(temp2);
    size(rho(:,50));
    %Computng the extra term (Gong ahead N_i steps)

    for j = i+1 : i+N_i-1
        mul_ext = x(:,j)' * Q * x(:,j) +  ( x(:,j)' * (K'*R*K) * x(:,j)) ;
    end
    
    %eta 50
    eta(1,50) = ( x(:,i)' * Q * x(:,i) ) + ( (K*x(:,i))' * R * (K*x(:,i)) ) + mul_ext - (x(:,i+N_i)' * (G12 * inv(G22) * G12') * x(:,i+N_i) ) ;
    
    if mod(i,50) == 0
            
            Z = rho';
            m = Z' * Z ;
            q = Z' * eta' ;
            o = rank(m);
            vG = inv(m) * q;
            
            f = size(vG);
            G11(1,:) = vG(1:4);
            G11(2,:) = vG(5:8);
            G11(3,:) = vG(9:12);
            G11(4,:) = vG(13:16);
            G11;
            G12(1,:) = vG(17:18);
            G12(2,:) = vG(19:20);
            G12(3,:) = vG(21:22);
            G12(4,:) = vG(23:24);
            G12;
            G22(1,:) = vG(25:26);
            G22(2,:) = vG(27:28);
            G22;
            K =  -inv(G22)*G12';
            i;
    end

    


end

G_z = [G11 G12;
       G12' G22];
G_z
figure(1);
sgtitle('State trajectories')
hold on;
subplot(2,2,1);
plot(x(1,:),'linewidth',1.2);
xlabel('Time step');
ylabel('pitch angle theta ');
hold on; 
subplot(2,2,2)
plot(x(2,:))
xlabel('Time step');
ylabel('Yaw angle ');
hold on; 
subplot(2,2,3);
plot(x(3,:))
xlabel('Time step');
ylabel('theta dot (rad)');
hold on; 
subplot(2,2,4);
plot(x(4,:))
xlabel('Time step');
ylabel('Yaw angle dot');
hold on;  
figure





