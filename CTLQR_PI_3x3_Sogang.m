%% Information
% Code from http://www.uta.edu/utari/acs/code/Software%20from%20Research.htm
% Related Article : D. Vrabie, O. Pastravanu, M. Abu-Khalaf, and F. L. Lewis,
% ¡°Adaptive optimal control for continuous-time linear systems based on policy iteration,¡±
%  Automatica, vol. 45, pp. 477-484, 2009.
%  Editted by MS Chang, Ph.D candidate from South Korea (19. 6. 3)

clear all;close all;clc;

%% Initialization
global P u T Tfinal epsilon A B C R HH  C1 C2
global zero_G pole_G gain_G ke kd kp wn ode_k noise

% sample time
T=0.05; Tfinal = 6;
N=Tfinal/T; % Length of the simulation in samples
epsilon = 0.00000001; % For the converge
R = 0.7;
HH = 10; % Period of update policy

%system matrix
zero_G = []; pole_G = [-5]; gain_G = 5; % Human's model (zero, pole, gain)
G_fh = zpk(zero_G,pole_G,gain_G); [num, den] = tfdata(G_fh, 'v'); 
ke = num(2); kp = den(2); kd = den(1); % Human's parameter

wn = 1; % Natural frequency
A = [0         1          0;       % System matrix
     -2*wn     -wn^2       0;
%      0     0       0;
     ke/kd     0          -kp/kd]; 
B = -[0; 1; 0];                     % Input matrix

C = [1 0 0];

eig_A = eig(A)


%initial conditions
x0_set = [10,0,10,0]
x0=x0_set; % 3 vectors and J
P=[0 0 0;
   0 0 0; 
   0 0 0]; % Positive definite and symmetric matrix P
Pold = eye(3);  Psave = zeros(length(P),length(P),2);  % For monitoring P matrix
uu=[];       % saving the control signal for plot
xx=[];       % saving the state for plot
KK=[];       % saving the optimal gain for plot

% Vectorization: Parameters returned by the least squares
WW=[P(1,1); 2*P(1,2); 2*P(1,3); P(2,2); 2*P(2,3); P(3,3)];
WWP=[WW; 0];

% Parameters for the batch least squares
j=0; Xpi=[];
E=real(eig(A-B*inv(R)*B'*P)); % saves the poles of the closed loop system
EE = [E];
upd=[];                   % stores information relative to updates the parameters
k=1; ch=0;                 % Real time iteration, Interval iteration
qm = zeros(10,1);
qd = zeros(10,1); qd(1) = x0_set(1);
qdd = zeros(10,1);
noise_save = 0;

%% Solving for the cost using least squares
figure('color','w');
while(1)
    
    j=j+1; % Iteration    
%     noise = rand(1)*0.5-0.25;
    noise = 0;
    noise_save(k) = noise;
    X(j,:)=[x0(1)^2 x0(1)*x0(2) x0(1)*x0(3)...
            x0(2)^2 x0(2)*x0(3)...
            x0(3)^2]; % Previous state vectors
    before_cost=[x0(1) x0(2) x0(3)]*P*[x0(1) x0(2) x0(3)]'; % Previous cost fn
    
    % Ode function simulation (Continuous time)
    ode_k = 0;
    tspan=[0 T];     [t,x]= ode45(@odefile,tspan,x0);

    % Update state vectors
%     if k>10, x(:,3) = x(:,3)-x(:,3)-2; end
    x1=x(length(x),1);     x2=x(length(x),2);     x3=x(length(x),3);
    xx = [xx; x];
    
    % Desired output
    if k>1,     qd(k) = 3*sin(0.7*pi*T*(k-1))+x0_set(1); qdd(k) =(qd(k)-qd(k-1))/T; end
    qm(k) = qd(k) - C*x(length(x),1:3)';
    
    % Prepare for the least square problem
    Xpi(j,:)=X(j,:)-[x1^2 x1*x2 x1*x3 x2^2 x2*x3 x3^2];
    Y(j,:)=x(length(x),4);
    x0=[x(length(t),1) x(length(t),2) x(length(t),3) 0];
    
    % Update cost fn
    after_cost=x(length(x),4)+[x1 x2 x3]*P*[x1;x2;x3];
    
    % Wait for storing the data
    if (abs(after_cost-before_cost)>0.00001)&&(ch==0)&&(mod(j,HH)~=1),
        j=0;
        ch=ch+1;
    else
        if abs(after_cost-before_cost)>0.00001,
            ch=ch+1;
        end
    end
    
    % the batch least squares is made on 10 values
    if mod(j,HH)==0,
        
        if (abs(after_cost-before_cost)>0.00001)&&(ch==HH),
            weights=(Xpi'*Xpi)\(Xpi'*Y); %calculation of the weights
%             JJ = Xpi'*Xpi;
%             [JJ_U JJ_D JJ_V]=svd(JJ);
%             JJ_DT = diag(diag(1./JJ_D));
%             JJinv = JJ_U*JJ_DT*JJ_V';
%             weights = JJ_U*JJ_DT*JJ_V'*(Xpi'*Y);
            
            upd=[upd 1];
            KK = [KK; -inv(R)*B'*P;];
        else
            %there is no reason to update
            upd=[upd 0];
            KK = [KK; -inv(R)*B'*P;];
        end
        
        WWP=[WWP [weights; k*T]];  WW=[WW weights];
        
        % Calculating P
        P=[weights(1)   weights(2)/2  weights(3)/2;
            weights(2)/2 weights(4)   weights(5)/2;
            weights(3)/2 weights(5)/2 weights(6)   ];
        
        % Stop the PI
%         if (norm(P-Pold) < epsilon && x1 < 0.01) , break; end % Stop the PI
        E = eig(A-B*inv(R)*B'*P);
        j=0;  ch=0;  EE=[EE real(E)];
        Pold = P; % Save the P
        Psave(:,:,k) = P; % For Monitoring P
        subplot(212);
        plot(0.05+T*(k-1),Psave(1,1,k),'b.','linewidth',2); hold on;
        plot(0.05+T*(k-1),Psave(1,2,k),'g+','linewidth',2); hold on;
        plot(0.05+T*(k-1),Psave(2,3,k),'ro','linewidth',2); hold on;
        plot(0.05+T*(k-1),Psave(3,3,k),'cs','linewidth',2); hold on;
    end % Finishi the if-then
    if k>N, break; end % Stop the PI
    
    subplot(211);
    plot(t+T*(k-1),x(:,1),'b','linewidth',2); hold on;
    plot(t+T*(k-1),x(:,2),'g:','linewidth',2);  hold on;
    plot(t+T*(k-1),x(:,3),'r-.','linewidth',2); hold on;
    
    subplot(212); hold on;
    drawnow;
    
    uu=[uu u]; % Monitoring for u
    k= k+1;
end % End of while loop
figure(1); 
subplot(211); legend('1','2','3')
% % axis([0 Tfinal -0.2 0.2]); title('System states'); ylabel('x(t)')
% % subplot(212); axis([0 Tfinal 0 max(max(P))*1.1]); set(gca,'box','on'); ylabel('P')
% xlabel('Time(s)');

%% Figure
t_fig = 0:T:T*(length(uu)-1);
sol=care(A+rand(size(A)),B,C'*C,1); % Solution of LQR by Dynamic programming 

% Controller output
figure('color','w');  
subplot(311);
plot(t_fig,uu-noise_save(1:end-1),'b.','linewidth',2); hold on;
% plot(t_fig,noise_save(1:end-1),'r.','linewidth',2); hold on;
ylabel('u {Control Input}');
if max(uu) == 0, uu(end) = 1; end
axis([0 Tfinal min(uu) max(uu)]); 
subplot(312);
plot(WWP(7,1:end-1),upd,'bo','linewidth',2);
ylabel('Policy update'); xlabel('Time(s)');
% axis([0 Tfinal -0.5 1.5])
subplot(313);
plot(WWP(7,1:end-1),KK(:,1),'b','linewidth',2); hold on;
plot(WWP(7,1:end-1),KK(:,2),'g','linewidth',2); hold on;
plot(WWP(7,1:end-1),KK(:,3),'r','linewidth',2); hold on;
plot(WWP(7,1:end-1),KK(:,1),'bo','linewidth',2); hold on;
plot(WWP(7,1:end-1),KK(:,2),'go','linewidth',2); hold on;
plot(WWP(7,1:end-1),KK(:,3),'ro','linewidth',2); hold on;
ylabel('Optimal gain K'); xlabel('Time(s)');
legend('K(1): K','K(2): B','K(3): kh')
% axis([0 Tfinal -0.5 1.5])
drawnow;

% P matrix parameter
figure('color','w');
plot(WWP(7,:),WWP(1:6,:)','linewidth',2); hold on;
plot(WWP(7,:),WWP(1:6,:)','ko','linewidth',1);

title('P matrix parameters'); xlabel('Time(s)')
% axis([0 Tfinal -0.5 5.5])
drawnow;

% Pole assignment
figure('color','w');
for jj=1:N/HH+1,
    plot(HH*T*(jj-1),EE(1,jj),'bo'); hold on;
    plot(HH*T*(jj-1),EE(2,jj),'r*'); hold on;
    plot(HH*T*(jj-1),EE(3,jj),'ks'); hold on;
end
legend('1','2','3')
title('Poles of the closed loop system'); xlabel('Time(s)')



% Robot closed loop System Modeling
numr = [719.3  1.208e04]; denr = [1 101.2 1956 4805];
Gs = tf(numr,denr);
Gz = c2d(Gs,T,'zoh'); [a b]=tfdata(Gz,'v')
Kp=5; Kd=0.005; C=tf([Kd Kp],[1]); Cz=c2d(C,T,'matched')
Gcz=feedback(Gz*Cz,1);
pole(Gcz)

Nr = length(qm);
y = zeros(1,Nr);
e = zeros(1,Nr);
ur = zeros(1,Nr);
urp = zeros(1,Nr);
urd = zeros(1,Nr);
tr = T*Nr;
y(1:3) = qm(1:3);
for k = 4:N,
    % Simulation by G(z), i.e., actual dynamics
    y(k) = -b(2)*y(k-1)-b(3)*y(k-2)-b(4)*y(k-3)+a(1)*ur(k)+a(2)*ur(k-1)+a(3)*ur(k-2)+a(4)*ur(k-3);
    e(k) = qm(k)-y(k);

    % Feedback control
    urp(k) = Kp*e(k);   % PD control
    urd(k) = Kd*(e(k)-e(k-1))/T;   % PD control
    ur(k) = urp(k) + urd(k);   % PD control
    
    % Saturation of control input
    if     ur(k) > 10,       ur(k) = 10;
    elseif ur(k) < -10,      ur(k) = -10;
    end
end

% Trajectory
figure('color','w');
subplot(311);
plot(t_fig,qd(1:end-1),'b','linewidth',2); hold on;
plot(t_fig,qm(1:end-1),'r:','linewidth',2); hold on;
plot(t_fig,y(1:end-1),'k--','linewidth',1); hold on;
ylabel('y'); xlabel('time (sec)')
legend('q_d','q_m','q')

subplot(312);
plot(t_fig,e(1:end-1),'b','linewidth',2); hold on;
plot(t_fig,qd(1:end-1)-qm(1:end-1),'r:','linewidth',2); hold on;
ylabel('err'); xlabel('time (sec)')

subplot(313);
plot(t_fig,ur(1:end-1),'b','linewidth',2); hold on;
plot(t_fig,ur(1:end-1),'r:','linewidth',2); hold on;
plot(t_fig,ur(1:end-1),'k--','linewidth',2); hold on;
ylabel('u (V)'); xlabel('time (sec)')
legend('u_p','u_d','u')


%% Ode45 (Realtime)
function xdot=odefile(t,x)
global P u A B C R ode_k C1 C2 noise;

Q = 1;

x1=x(1);
x2=x(2);
x3=x(3);
%calculating the control signal
% u=0;
u=-inv(R)*B'*P*[x1;x2;x3]+noise;
% if u > 10, u = 10; elseif u < -10, u = -10; end

%updating the derivative of the state=[x(1:2) V]
xdot=[A*[x1;x2;x3]+B*u
      [x1;x2;x3]'*Q*[x1;x2;x3]+u'*R*u];
  ode_k = ode_k+1;
end
