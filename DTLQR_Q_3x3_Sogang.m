close all; clc; clear all; warning off;
global dt Tfinal gamma R Q N n epsilon
%% Offline Initialization
epsilon = 0.7; % Most important for convergence!! (Learning Stop)
Tfinal = 0.5; dt = 0.001; 
t = 0:dt:Tfinal; Nt = length(t);
% Gs = zpk([],[2 5],1);
% Gz = c2d(Gs,dt,'zoh'); [numz denz]=tfdata(Gz,'v');
% [A B C D] = tf2ss(numz,denz);
% C = [10 0];

m = 1; c = 23; k = 50;
F = [0 1; -k/m -c/m]; Gb = [0; 1000000/m];
C = [1 0]; D =0;
Gc = ss(F,Gb,C,D);
Gd = c2d(Gc,dt,'zoh');
[A, B, C, D] = ssdata(Gd);

F=[-0.9]; L=[0 0 0];
R = 0.3;  Q = 1;
gamma = 0.8;
A1=[A    [0;0];
    [0 0] F];
B1=[B;0];
% C1 = [C 0; zeros(1,length(C)+1); zeros(1,length(C)+1); ]-eye(length(C)+1);
% Q1=C1'*Q*C1;
Q1=[C'*Q*C -C'*Q;-Q*C Q];
G=[Q1 zeros(length(Q1),1);
    zeros(1,length(Q1)) R];

x=[10;0]; r=[10]; % Initial value
X=[x; r];
X_off=[x; r];

%% Offline Optimization
while(1)
    % offline solution
    P1=dare(sqrt(gamma)*A1,sqrt(gamma)*B1,Q1,R);   % P by LQR
    H1 = [Q1 + gamma*A1'*P1*A1   gamma*A1'*P1*B1;
        gamma*B1'*P1*A1        R+gamma*B1'*P1*B1];
    
    H1yy=H1(length(H1),length(H1)); H1yx=H1(length(H1),1:length(H1)-1);
    K1=-inv(H1yy)*H1yx;
    
    % offline Simulation
    for k=1:10
        H1=G+gamma*[A1 B1;K1*A1 K1*B1]'*H1*[A1 B1;K1*A1 K1*B1];
        H1yy=H1(length(H1),length(H1)); H1yx=H1(length(H1),1:length(H1)-1);
        K1=-inv(H1yy)*H1yx;
    end
    
    for k=1:Nt
        u_off(k)=K1*X_off(:,k);
        X_off(:,k+1)=A1*X_off(:,k)+B1*u_off(k);
        y_off(k)=C*X_off(1:2,k);
    end
    e_off = X_off(3,1:end-1)-y_off;
    figure('color','w')
    subplot(311);
    plot(t,X_off(3,1:end-1),'b','linewidth',2); hold on;
    plot(t,y_off,'r:','linewidth',2); hold on;
    legend('r','y'); ylabel('Position');
    subplot(312);
    plot(t,e_off,'b','linewidth',2);
    ylabel('e'); xlabel('Time (s)')
    subplot(313);
    plot(t,u_off,'b','linewidth',2);
    ylabel('u'); xlabel('Time (s)')
    drawnow;
    break;
end


%% online solution
N = 20;
H = eye(length(A1)+1);
Hold = 10*eye(length(A1)+1);
Hyy=H1(4,4);  Hyx=H1(4,1:3);  K=-inv(Hyy)*Hyx;
n = length(H)*(length(H)+1)/2;
% H = H1;
H = eye(4)*0.01;
Hyy=H(4,4); Hyx=H(4,1:3); K=-inv(Hyy)*Hyx;
% K = [-1 0 0]


zbar = zeros(n,N); d_target = zeros(N,1); kk = 1;
Ysave = []; Xsave =[]; Ksave = [];
i=1; isave =0; Xpi = []; Y = []; Z = []; d=[]; d1=[]; d2=[];
j=1; h=0; ranksave = []; update = 0; detsave = []; noise = 0;
u = zeros(Nt,1);
figure('color','w');

%% Iteration Start !!

% Reference
for k=1:2*Nt+1
%         r(k) = 1*cos(50*pi*dt*(k-1));
end

while(1)
    X(:,i)=[x(:,i);r(:,i)]; % Current State
    r(:,i+1)=F*r(:,i);      % Next Trajectory
    
    % Policy Update (Tricky for nonsingularity)
    noise=0.01; if t(i)> 0.5, noise = 0; end
    BB(i) = noise*rand(1); % Adding noise to avoid singularity (u is dependent on x)
    u(i)=K*X(:,i) + BB(i); % Noisy input
    Z(:,i)=[X(:,i); u(i)]; % State for Q function approximation 
    
    % System model
    x(:,i+1)=A*x(:,i)+B*u(i);
    y(i)=C*x(:,i);
    X(:,i+1)=[x(:,i+1);r(:,i+1)];  % New State for Q function
    
    % Target
    d_target=[X(:,i); u(i)]'*G*[X(:,i); u(i)]+gamma*[X(:,i+1);K*X(:,i+1)]'*H*[X(:,i+1);K*X(:,i+1)];
    zbar=[X(1,i)^2; X(1,i)*X(2,i); X(1,i)*X(3,i); X(1,i)*u(i); X(2,i)^2; X(2,i)*X(3,i); X(2,i)*u(i); X(3,i)^2; X(3,i)*u(i); u(i)^2];
    
    % Initialization for the Least Square
    if h == 1 && i<=Tfinal/dt,
        h = 0;         Xpi = zeros(n,N);        Y = zeros(N,1);
    end
    
    % Collect target during N steps
    Xpi(:,i-isave) = zbar;     Y(i-isave,:) = d_target;
    Xsave(i,:) = zbar;         Ysave(i,:) = d_target;
    if i-isave > 3,  Tsave(i,:) = [mod(i,N), Xpi(1,3)];    end
    
    % Learning & Least square problem
    if mod(i,N) == 0
        if i>Tfinal/dt, break; end
        h = 1;        eL=abs(K-K1);        kk = kk+1;
        update = 0;
        if norm(H-Hold) > epsilon
%            epsilon = norm(H-Hold)*10;
            %             if i>Tfinal/dt || rank(Xpi*Xpi') ~= n, sprintf('Rank(%.4f) Error, Det = 0;',rank(Xpi*Xpi')), break; end
            ranksave = [ranksave rank(Xpi*Xpi')]; % Check the rank for the singularity 
            detsave = [detsave det(Xpi*Xpi')];
            vH=(Xpi*Xpi')\(Xpi*Y);                % New vectorization of H
            Hold = H; 
            H=[vH(1,1) vH(2,1)/2 vH(3,1)/2 vH(4,1)/2 ;  % New H
                vH(2,1)/2 vH(5,1) vH(6,1)/2 vH(7,1)/2;
                vH(3,1)/2 vH(6,1)/2 vH(8,1) vH(9,1)/2;
                vH(4,1)/2 vH(7,1)/2 vH(9,1)/2 vH(10,1)];
            
            Hyy=H(4,4);Hyx=H(4,1:3);
            K=-inv(Hyy)*Hyx;
            update = 1;
        end
        d(kk) = norm(eL);         d1(kk) = norm(H1-H);
        d2(kk) = norm(H-Hold);    d3(kk) = update;
        Ksave(kk,:) = K; % Save the New optimal gain
        
        j=j+1;
        isave = i;
    end
    
    % Realtime Monitoring
    %     if i>2
    %         plot(t(i-1:i),r(i-1:i),'b','linewidth',2); hold on;
    %         plot(t(i-1:i),y(i-1:i),'r:','linewidth',2);
    %         drawnow;
    %     end
    
    if i>Tfinal/dt , sprintf('Success'), break; end
    i=i+1;
end

% For Figure
rf = r(1:Nt); yf = y(1:Nt); uf = u(1:Nt); ef = rf-yf; t1= linspace(0,Tfinal,kk);

% Overall Monitoring
figure('color','w');
subplot(211);
plot(t,rf,'b','linewidth',2); hold on;
plot(t,yf,'r:','linewidth',2);
plot(t,y_off,'g--','linewidth',2); hold on; % Offline LQR
legend('r','y_{Online}','y_{Offline}');
ylabel('Output'); xlabel('Time(s)')
% subplot(312);
% plot(t,ef,'b','linewidth',2);
% ylabel('e'); xlabel('Time (s)')
subplot(212);
plot(t,uf,'b.','linewidth',2);
ylabel('u'); xlabel('Time (s)')

% Update Monitoring
figure('color','w');
% subplot(311);
plot(t1,d3,'bo','linewidth',2); hold on;
plot(t1,d3,'r','linewidth',2);
ylabel('Policy Update'); xlabel('Time (s)')
% subplot(212);
% plot(t1,Ksave(:,1),'b','linewidth',2); hold on;
% plot(t1,Ksave(:,2),'r','linewidth',2); hold on;
% plot(t1,Ksave(:,3),'k','linewidth',2); hold on;
% legend('K(1)','K(2)','K(3)','','','')
% plot(t1,Ksave(:,1),'bo','linewidth',2); hold on;
% plot(t1,Ksave(:,2),'ro','linewidth',2); hold on;
% plot(t1,Ksave(:,3),'ko','linewidth',2); hold on;
% ylabel('K'); xlabel('Time (s)')

% Value and policy monitoring
figure('color','w');
subplot(211);
plot(t1,d,'bo','linewidth',2); hold on
plot(t1,d,'k','linewidth',2);
ylabel('|| K_{RL}-K_{LQ} ||')
% subplot(312);
% plot(t1,d1,'bo','linewidth',2); hold on
% plot(t1,d1,'r','linewidth',2)
% ylabel('|| H_{RL}-H_{LQ} ||')
subplot(212);
plot(t1,d2,'bo','linewidth',2); hold on
plot(t1,d2,'k','linewidth',2)
ylabel('|| H_{j+1}-H_{j} ||')
xlabel('Time (s)')
