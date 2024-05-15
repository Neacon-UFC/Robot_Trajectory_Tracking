clear all
close all
clc


Ur = [[0.3;0.15],[0.3;-0.15]];
ubar = Ur;
K3 = 0;
K4 = 0;
stop2 = 0;

tic
for j=1:length(Ur)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       SUM-OF-SQUARES FEASIBILITY (OR OPTIMIZATION) PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initializations
pvar e1 e2 e4 e5
e_n=[e1;e2;e4;e5];
ubar = Ur;

% Polynomial recasted robot tracking model from eq.(16).
f = [Ur(1,j)*(e5+1); Ur(1,j)*e4; Ur(2,j)*(e5+1); -Ur(2,j)*e4];
g = [-1 e2; 0 -e1; 0 -e5-1; 0 e4];

% Initial V to be fixed in the first program of the iteration 1.
Vo = e_n'*eye(length(e_n))*e_n;
V = Vo;

% Maximum iteration.
i = 0;
i_max = 100;

% Monomials of the controllers [v;ω] whose each coefficient (the controller
% gains) will be found as part of the SOS Program solution.
z = monomials(e_n,1:4);

% Tolerance β between the previous and current solution
beta = 0.1;

% Inicial choice for the controller coefficients to enable the first check
% of the stopping criterion. We do set a large value to all coefficient so the first check
% be guaranteed to not fall within the tolerance β.
uo = ones(length(monomials(e_n,1:4)),size(g,2))*100;

K11 = K3;
K22 = K4;
stop1 = stop2;

% Factibility loop until solution convergence.
while i<i_max
    %%%%% SOS-PROGRAM 1: SET OF CONTROLLER GAINS SOLUTION FROM THE PREVIOUS FIXED V %%%%%
    % solving with SOSTOOLS (V fixed).
    % Defining the first set of decision variables (u,l2,s1).
    prog = sosprogram(e_n);
    [prog,K1] = sospolymatrixvar(prog,monomials(e_n,0),[1,length(z)]);
    [prog,K2] = sospolymatrixvar(prog,monomials(e_n,0),[1,length(z)]);
    u=[K1*z;K2*z];
    [prog,l2] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)],'wscoeff');
    % Defining the SOS constraints
    g1 = e4^2+e5^2-1; % constraint from eq.(13)
    gradV = [diff(V,e1);diff(V,e2);diff(V,e4);diff(V,e5)];
    dV = -gradV'*(f+g*(u+ubar(:,j)))-s1*g1-l2;
    prog = sosineq(prog,dV);
    % First solution: (u,l2,s1)
    sol = sossolve(prog);
    u = sosgetsol(sol,u);
    K1 = double(sosgetsol(sol,K1));
    K2 = double(sosgetsol(sol,K2));
    l2 = sosgetsol(sol,l2);
    s1 = sosgetsol(sol,s1);
    %%%%% SOS-PROGRAM 2: LYAPUNOV FUNCTION SOLUTION FROM THE PREVIOUS CONTROLLER GAINS %%%%%
    % solving with SOSTOOLS (u fixed)
    % Defining the second set of decision variables (V,,l1,l2,s1).
    prog = sosprogram(e_n);
    [prog,V] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,l1] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,l2] = sospolyvar(prog,[monomials(e_n,2:4)],'wscoeff');
    [prog,s1] = sospolyvar(prog, [monomials(e_n,2:4)], 'wscoeff');
    % Defining the SOS constraints
    g1 = e4^2+e5^2-1; % constraint from eq.(13)
    gradV = [diff(V,e1);diff(V,e2);diff(V,e4);diff(V,e5)];
    dV = -gradV'*(f+g*(u+ubar(:,j)))-s1*g1-l2;
    prog = sosineq(prog,dV);
    prog = sosineq(prog,V-l1);
    % Second solution: (V,l1,l2,s1)
    sol = sossolve(prog);
    V = sosgetsol(sol,V);
    l1 = sosgetsol(sol,l1);
    l2 = sosgetsol(sol,l2);
    s1 = sosgetsol(sol,s1);

    % Stopping criterion
    delta = max(abs(full(u.coefficient)-uo))
    if delta<beta
        break
    end
    % Stores the current controller coefficients to be compared in the next
    % loop.
    uo = full(u.coefficient);
    i = i+1
end

stop2 = i;
K3 = K1;
K4 = K2;

end
toc

%delta_V = max(abs(full(V.coefficient)-Vo)) in case
%Vo = full(V.coefficient);


%%
% The next loop difines the reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Geracao da trajetoria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.1;          % sampling time
h = 0.1;           % integrating time

xr = [0.7;1.5;-pi/8];  % x,y,theta referenciais para in�cio de trajet�ria. 
Xr = xr;         %Euler

Tsim = 42;

for k = 1:round(Tsim/h)  % 1 a tempo de simula��o/ tempo de integra��o   
    vr(k) = 0.3;
    if k<=round(Tsim/(2*h))
        wr(k) = 0.15;
     else
         wr(k) = -0.15;
    end
    %Xr0 = [Xr0 xr0+h*([vr(k)*cos(xr0(3));vr(k)*sin(xr0(3));wr(k)])];  %Euler
    %Xr = [Xr xr+h*([vr(k)*sin(xr(3));vr(k)*cos(xr(3));wr(k)])];  %Euler
    Xr = [Xr xr+h*([vr(k)*cos(xr(3));vr(k)*sin(xr(3));wr(k)])];  %Euler
    xr = Xr(:,k+1);
end
plot(Xr(1,1),Xr(2,1),'o','MarkerFaceColor',	"#4DBEEE");
hold on
plot(Xr(1,:),Xr(2,:),'LineWidth',2,'Color',	"#4DBEEE");
hold on

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          ROBOT REFERENCE TRACKING SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%xo = [0.8;-0.3;pi/4];
xo = [0;1;0];
x = xo;
X = x;
v = [];
w = [];

syms E1 E2 E4 E5
E = [E1;E2;E4;E5-1];
Z = monomials(E,1:4);

for k = 1:(round(Tsim/h))    % 1 a tempo de simulacao/tempo de integracao  
    
    a = Xr(:,k)-x;
    E1 = cos(x(3))*a(1)+sin(x(3))*a(1);
    E2 = -sin(x(3))*a(1)+cos(x(3))*a(2);
    E3 = a(3);
    E4 = sin(E3);
    E5 = cos(E3);
    
    if k<=round(Tsim/(2*h))
        v(k) = K11*subs(Z)+ubar(1,1);
        w(k) = K22*subs(Z)+ubar(2,1);
    else
        v(k) = K3*subs(Z)+ubar(1,2);
        w(k) = K4*subs(Z)+ubar(2,2);
    end

    X = [X x+dt*[v(k)*cos(x(3));v(k)*sin(x(3));w(k)]];
    x = X(:,k+1);
end

Error = Xr(:,1:end-1)-X(:,1:end-1);
Error4 = sin(Error(3,:));
Error5 = cos(Error(3,:));



plot(X(1,1),X(2,1),'o','MarkerFaceColor',"#FF0000")
hold on
plot(X(1,1:end-1),X(2,1:end-1),'Color',"#FF0000")
legend('q_{r}(0)','Reference Robot','q_{c}(0)','Real Robot')
xlabel('x (m)')
ylabel('y (m)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotagem de v e w 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
time = dt*(0:Tsim/h-1);
subplot(2,1,1)
plot(time,v,'b','LineWidth',1.5)
xlabel('Time (s)')
ylabel('U1 (m/s)')
legend('v')

hold on
subplot(2,1,2)
plot (time,w,'r','LineWidth',1.5)
xlabel('Time (s)')
ylabel('U2 (rad/s)')
legend('ω')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Plotagem do erro x y e theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(time,Error(1,:),time,Error(2,:),time,Error(3,:),time,Error4(:),time,Error5(:),'LineWidth',1.5);
legend('e1','e2','e3','e4','e5')
