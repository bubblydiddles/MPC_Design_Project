%%
clear
close
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%% Defining the robot dynamics
R = [-sqrt(3)/2 0.5 1;
    0 -1 1;
    sqrt(3)/2 0.5 1];
A = zeros(3);
B = inv(R);
C = eye(3);
D = 0;

Ts = 0.1; %Sampling Time
T = 10; %Total Simulation Time
t_span = 0:Ts:T;
x0 = [2, 2 , pi/2]';
simulation_steps = T/Ts; 

robot = ss(A,B,C,D); %Continuous Dynamics
discrete_robot = c2d(robot, Ts); %Discretized Dynamics

u_bound = [-1 1];
x_upper_bound = [10, 10, 2*pi]';
x_lower_bound = [-10, -10, 0]';

%% Infinite Horizon LQR Control
Q = 0.1*eye(3);
R = 0.01*eye(3);

[lqr_controller,~,~] = dlqr(discrete_robot.A, discrete_robot.B, Q, R);
lqr_controlled_robot = ss((robot.A-robot.B*lqr_controller), zeros(3), C, D);

x = lsim(lqr_controlled_robot, zeros(length(t_span),3), t_span, x0)';

PLOTTER(x);

%% Plotting the Terminal Set
bounds = ones(6,1);
plotregion([-lqr_controller; lqr_controller],-bounds)

%% MPC Problem Definition
dim.nx = 3; %State Dimension
dim.nu = 3; %Input Dimension
dim.N = 50; %Prediction Horizon

weight.Q = 0.1*eye(3);
weight.R = 0.01*eye(3);

weight.P = dare(discrete_robot.A, discrete_robot.B, weight.Q, weight.R);

prediction_model = predmodgen(discrete_robot, dim);
[Hx,hx] = costgen(prediction_model, weight, dim);

Hx = (Hx+Hx')/2;

%% Point-To-Point Controller(State Feedback, No disturbances)
x = zeros(dim.nx,simulation_steps+1);
x(:,1) = x0;

u1 = zeros(dim.nu, simulation_steps);

for i = 1:simulation_steps
    u_uncon = sdpvar(dim.nu*dim.N,1);
    bound = ones(length(u_uncon), 1);
    Constraint = [u_bound(1)*bound <= u_uncon<= u_bound(2)*bound, kron(ones(dim.N+1,1), x_lower_bound) <= prediction_model.T*x(:,i) + prediction_model.S*u_uncon <= kron(ones(dim.N+1,1), x_upper_bound), u_bound(1)*ones(3,1)  <=lqr_controller*(prediction_model.T(end-2:end,:)*x(:,i) + prediction_model.S(end-2:end,:)*u_uncon)<= u_bound(2)*ones(3,1)];
    Objective = 0.5*u_uncon'*Hx*u_uncon + (hx*x(:,i))'*u_uncon;
    optimize(Constraint, Objective);
    u_uncon = value(u_uncon);
    
    u1(:,i) = u_uncon(1:dim.nu)';
    x(:,i+1) = x(:,i) + robot.B*Ts*u1(:,i);
end

PLOTTER(x);

%% Observer Design

Bd = discrete_robot.B*[1;0;0];
Cd = [0; 0; 0];
Atilde = [discrete_robot.A Bd;0 0 0 1];
Btilde = [discrete_robot.B ; zeros(1,3)];
Ctilde = [discrete_robot.C Cd];

A_dist = [eye(3)-discrete_robot.A -Bd;discrete_robot.C Cd];
%% Disturbed P2P MPC definition
d=0.2;
dime.nx = 4; %State Dimension
dime.nu = 3; %Input Dimension
dime.nd = 1; %Disturbance Dimension
dime.N = 50; %Prediction Horizon

weighte = weight;
weighte.Q = blkdiag(weighte.Q, zeros(dime.nd));
weighte.P = blkdiag(weighte.P, zeros(dime.nd));

LTIE.A = Atilde;
LTIE.B = Btilde;
LTIE.C = Ctilde;

prediction_model = predmodgen(LTIE, dime);
[Hx,hx] = costgen(prediction_model, weighte, dime);

Hx = (Hx+Hx')/2;

%% P2P Offset-Free MPC(Output-Feedback)

poles = [0.8 , 0.9 , 0.6, 0.5];
L=(place(LTIE.A',LTIE.C',poles))';

x = zeros(dime.nx,simulation_steps+1);
x(:,1) = [x0; d];

xhat = zeros(dime.nx,simulation_steps+1);
xhat(:,1) = x(:,1) + [0.1, 0, 0, 0]';

for i = 1:simulation_steps
    u_uncon = sdpvar(dime.nu*dime.N,1);
    bound = ones(length(u_uncon), 1);
    Constraint = [u_bound(1)*bound <= u_uncon<= u_bound(2)*bound,  kron(ones(dime.N+1,1), [x_lower_bound; -Inf]) <= prediction_model.T*x(:,i) + prediction_model.S*u_uncon <= kron(ones(dime.N+1,1), [x_upper_bound; Inf]), u_bound(1)*ones(3,1)  <=lqr_controller*(prediction_model.T(end-3:end-1,:)*x(:,i) + prediction_model.S(end-3:end-1,:)*u_uncon)<= u_bound(2)*ones(3,1)];
    Objective = 0.5*u_uncon'*Hx*u_uncon + (hx*xhat(:,i))'*u_uncon;
    optimize(Constraint, Objective);
    u_uncon = value(u_uncon);
    
    u(:,i) = u_uncon(1:dime.nu)';
    x(:,i+1) = LTIE.A*x(:,i) + LTIE.B*u(:,i);
    
    xhat(:,i+1) = LTIE.A*xhat(:,i) + LTIE.B*u(:,i) + L*(x(1:3,i) - LTIE.C*xhat(:,i));
end

PLOTTER(x, xhat);
figure, hold on, plot(t_span(1:end-1), u(1,:)), plot(t_span(1:end-1), u(2,:)), plot(t_span(1:end-1), u(3,:)), legend('$\omega_1$', '$\omega_2$', '$\omega_3$'), title('Control Inputs');

%% Trajectory Generation

T = 5;
simulation_steps = T/Ts;

t_span = 0:Ts:T;
ur = zeros(length(t_span),3);
ur(:,1) = 1*sin(t_span);
ur(:,2) = 1*cos(2*t_span);
ur(:,3) = 1;

xr = lsim(robot, ur, t_span)';
xr(3,:) = mod(xr(3,:), 2*pi);

%% Trajectory Tracking - MPC Problem Definition
dim.nx = 3; %State Dimension
dim.nu = 3; %Input Dimension
dim.N = 20; %Prediction Horizon

weight.Q = 0.1*eye(3);
weight.R = 0.01*eye(3);

error_dynamics = ss(zeros(3), -robot.B, eye(3), zeros(3))';
discrete_error_dynamics = c2d(error_dynamics, Ts);

weight.P = dare(discrete_error_dynamics.A, discrete_error_dynamics.B, weight.Q, weight.R);

% Set these for simulation
init_error = [0 0.5 -pi/8];

prediction_model_e = predmodgen(discrete_error_dynamics, dim);
[He,he] = costgen(prediction_model_e, weight, dim);

He = (He+He')/2;


prediction_model = predmodgen(discrete_robot, dim);
[Hx,hx] = costgen(prediction_model, weight, dim);

Hx = (Hx+Hx')/2;

%% Trajectory Tracking
tracking_error = zeros(dim.nx,simulation_steps+1);
tracking_error(:,1) = init_error;

ub = zeros(dim.nu, simulation_steps);

for i = 1:simulation_steps
    u_uncon = sdpvar(dim.nu*dim.N,1);
    if i+dim.N<=simulation_steps
        uf = ur(i:i+dim.N-1,:);
    else
        uf = ur(i:end-1,:);
    end
    temp = reshape(uf.',1,[])';
    uf = zeros(size(u_uncon));
    uf(1:size(temp,1)) = temp;
    
    bound = ones(length(u_uncon), 1);
    Constraint = [ u_bound(1)*bound <= u_uncon+uf <= u_bound(2)*bound, kron(ones(dim.N+1,1), x_lower_bound) <= prediction_model.T*(xr(:,i)-tracking_error(:,i)) + prediction_model.S*(u_uncon+uf) <= kron(ones(dim.N+1,1), x_upper_bound)];
    Objective = 0.5*u_uncon'*He*u_uncon + (he*tracking_error(:,i))'*u_uncon;
    optimize(Constraint, Objective);
    u_uncon = value(u_uncon);
    
    ub(:,i) = u_uncon(1:dim.nu)';
    tracking_error(:,i+1) = tracking_error(:,i) + discrete_error_dynamics.B*ub(:,i);
end

x = zeros(size(xr));
x = xr - tracking_error;

u = ub+ur(1:end-1,:)';

figure, hold on, plot(xr(1,:), xr(2,:)), plot(x(1,:), x(2,:)), legend('Reference', 'Robot Coordinates'), title('XY Coordinates');
figure, hold on, plot(t_span, xr(3,:)),  plot(t_span, x(3,:)), legend('Reference', 'Robot Angle'), title('Angle');
figure, hold on, plot(t_span(1:end-1), u(1,:)), plot(t_span(1:end-1), u(2,:)), plot(t_span(1:end-1), u(3,:)), legend('$\omega_1$', '$\omega_2$', '$\omega_3$'), title('Control Inputs');
