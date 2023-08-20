close all;
clc;
clear all;

% System Paramerters (After Discretization)
A=[ 1      0.1   0.2955  0.01043
          0        1    5.073   0.2955
          0        0   0.2816  0.07465
          0        0   -12.33   0.2816]

B=[0.007831;
       0.1314;
     -0.02165;
      -0.3717]
 R=1;
Q=5*eye(4);

% Initial state
x=[0.1 0.1 0.1 0.1]';

G=[Q [0;0;0;0];[0 0 0 0] R];

% Arbitary Control Policy
K0=[0 0 0 0];
K=K0;

% Based on LQR (ARE)
P=idare(A,B,Q,R);

zbar(1:15,1:15)=0;
gamma(1:15,1)=0;
H_optimal=[Q+A'*P*A A'*P*B; B'*P*A  R+B'*P*B]

% Initailizing H to zero
H=[0 0 0 0 0 ;0 0 0 0 0 ;0 0  0 0 0;0 0 0 0 0;0 0 0 0 0 ];
H1=[];H2=[];H3=[];H4=[];H5=[];

for i=1:800
    a1=1;
 noise1= a1*0.01*(sin(1*i)^2*cos(9*i)+sin(2*i)^2*cos(0.1*i)+sin(-1.2*i)^2*cos(0.5*i)+sin(i)^5+sin(1.12*i)^2+cos(2.4*i)*sin(2.4*i)^3);
    
    % New Control Policy
    u=K*x(:,i)+noise1;  %probing noise added 
    
    % The new state
    x(:,i+1)=A*x(:,i)+B*u;

    % Collecting the data

    for j = 1:14
        gamma(j,1) = gamma(j+1,1);
    end


    % The new datapoint
    gamma(15,1)=[x(:,i); u]'*G*[x(:,i); u]+[x(:,i+1);K*x(:,i+1)]'*H*[x(:,i+1);K*x(:,i+1)];

    for j = 1:14
        zbar(:,j) = zbar(:,j+1);
    end

    temp = [x(1,i)^2; 
        x(1,i)*x(2,i); 
        x(1,i)*x(3,i) ;
        x(1,i)*x(4,i) ;
        x(1,i)*u;
        x(2,i)^2;
        x(2,i)*x(3,i);
        x(2,i)*x(4,i);
        x(2,i)*u;
        x(3,i)^2;
        x(3,i)*x(4,i);
        x(3,i)*u;
        x(4,i)^2;
        x(4,i)*u;
        u^2];
    
    % The new data point
    zbar(:,15)=temp;


    if mod(i,15)==0
       if i<=800
         m=zbar*zbar';  
         q=zbar*gamma;
         rank(m);
         vH=inv(m)*q;
         H=[vH(1,1) vH(2,1)/2 vH(3,1)/2 vH(4,1)/2 vH(5,1)/2; vH(2,1)/2 vH(6,1) vH(7,1)/2 vH(8,1)/2 vH(9,1)/2;vH(3,1)/2 vH(7,1)/2 vH(10,1) vH(11,1)/2 vH(12,1)/2;vH(4,1)/2 vH(8,1)/2 vH(11,1)/2 vH(13,1) vH(14,1)/2; vH(5,1)/2 vH(9,1)/2 vH(12,1)/2 vH(14,1)/2 vH(15,1)];                
         Huu=H(5,5);
        Hux=H(5,1:4);
       H5=[H5,H(5,5)];
        H1=[H1,H(5,1)];
        H2=[H2,H(5,2)];
        H3=[H3,H(5,3)];
        H4=[H4,H(5,4)];
        K=-inv(Huu)*Hux; 
       end
    end
 end

H

figure(1);
sgtitle('State trajectories')
hold on;
subplot(2,2,1);
plot(x(1,:),'linewidth',1.2);
xlabel('Time step');
ylabel('Linear Position x (m)');
hold on; 
subplot(2,2,2)
plot(x(2,:))
xlabel('Time step');
ylabel('Linear speed v (m/s)');
hold on; 
subplot(2,2,3);
plot(x(3,:))
xlabel('Time step');
ylabel('tilt angle 𝜃 (rad)');
hold on; 
subplot(2,2,4);
plot(x(4,:))
xlabel('Time step');
ylabel('tilt angular velocity 𝜔 (rad/s)');
hold on; 
figure

plot(H1,'-o','Color','b','LineWidth',1)
hold on;

plot(H2,'-square','Color','r','LineWidth',1)
hold on;

plot(H3,'-diamond','Color','y','LineWidth',1)
xlim([0 50])
ylim([-10 50])
hold on;

plot(H4,'-+','Color',[0.5 0 0.8],'LineWidth',1)
hold on;

plot(H5,'-X','Color','g','LineWidth',1)
xlabel('Iterations');
ylabel('Parameter Estimates');

title("Subsystem I")
hold on;

legend('Hux(1)','Hux(2)','Hux(3)','Hux(4)','Huu')




