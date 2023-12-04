% ------------------------------------------------------------------------
% Bending Planes Calibration
% ------------------------------------------------------------------------
% This MATLAB script reads trajectory data from two sets of points
% corresponding to antero-posterior (AP) and medio-lateral (ML) movements.
% It calculates the mean vectors for superior and inferior portions of the
% trajectory, estimates the initial plane normals, and then uses an
% optimization solver to refine these normals. Finally, it visualizes the
% trajectories and the optimized bending planes in a 3D plot.
%
% Author: Benjamin Fortuno

%%
close all
clear all

%% -----------------------------------------------------------------------
% ANTERO-POSTERIOR TRAJECTORY
% ------------------------------------------------------------------------
% Read data from an xlsx file for antero-posterior (AP) movement
Tap = readtable("motor2_current60.xlsx");

% Extract base and tip positions
bpx = Tap.BasePositionX';
bpy = Tap.BasePositionY';
bpz = Tap.BasePositionZ';
bs = [bpx; bpy; bpz];
mbs = mean(bs,2);

tpx = Tap.TipPositionX';
tpy = Tap.TipPositionY';
tpz = Tap.TipPositionZ';
tp = [tpx; tpy; tpz];
tb_vectors_ap = tp - mbs;

% Divide trajectory into superior and inferior portions
n = size(tb_vectors_ap,2);
allIndices = 1:n;
indicesSuperior = 1:n/2;
indicesInferior = setdiff(allIndices, 1:n/2);
tb_vectors_sup = tb_vectors_ap(:, indicesSuperior);
tb_vectors_inf = tb_vectors_ap(:, indicesInferior);
mtb_vectors_sup = mean(tb_vectors_sup, 2);
mtb_vectors_inf = mean(tb_vectors_inf, 2);

% Initial estimate of the normal to the AP bending plane
normal_ap0 = cross(mtb_vectors_sup, mtb_vectors_inf);
normal_ap0 = normal_ap0 / norm(normal_ap0);

% Visualization of the AP trajectory and initial normal
figure(1)
hold on
f(1) = plot3(tb_vectors_sup(1,:), tb_vectors_sup(2,:), tb_vectors_sup(3,:), 'Color','r','DisplayName', 'AP measured trajectory');
plot3([0 mtb_vectors_sup(1,:)], [0 mtb_vectors_sup(2,:)], [0 mtb_vectors_sup(3,:)], 'Color','r')
plot3(tb_vectors_inf(1,:), tb_vectors_inf(2,:), tb_vectors_inf(3,:), 'Color', 'r')
plot3([0 mtb_vectors_inf(1,:)], [0 mtb_vectors_inf(2,:)], [0 mtb_vectors_inf(3,:)], 'Color','r')
hold off

%% -----------------------------------------------------------------------
% MEDIO-LATERAL TRAJECTORY
% ------------------------------------------------------------------------
% Read data from an xlsx file for medio-lateral (ML) movement
Tml = readtable("motor1_current60.xlsx");

% Extract base and tip positions for ML movement
bpx = Tml.BasePositionX';
bpy = Tml.BasePositionY';
bpz = Tml.BasePositionZ';
bs = [bpx; bpy; bpz];
mbs = mean(bs,2);

tpx = Tml.TipPositionX';
tpy = Tml.TipPositionY';
tpz = Tml.TipPositionZ';
tp = [tpx; tpy; tpz];
tb_vectors_ml = tp - mbs;

% Divide ML trajectory into superior and inferior portions
n = size(tb_vectors_ml, 2);
allIndices = 1:n;
indicesSuperior = 1:n/2;
indicesInferior = setdiff(allIndices, 1:n/2);
tb_vectors_sup = tb_vectors_ml(:, indicesSuperior);
tb_vectors_inf = tb_vectors_ml(:, indicesInferior);
mtb_vectors_sup = mean(tb_vectors_sup, 2);
mtb_vectors_inf = mean(tb_vectors_inf, 2);

% Initial estimate of the normal to the ML bending plane
normal_ml0 = cross(mtb_vectors_sup, mtb_vectors_inf);
normal_ml0 = normal_ml0 / norm(normal_ml0);

% Visualization of the ML trajectory and initial normal
figure(1)
hold on
f(2) = plot3(tb_vectors_sup(1,:), tb_vectors_sup(2,:), tb_vectors_sup(3,:), 'Color','b', 'DisplayName', 'ML measured trajectory')
plot3([0 mtb_vectors_sup(1,:)], [0 mtb_vectors_sup(2,:)], [0 mtb_vectors_sup(3,:)], 'Color','b')
plot3(tb_vectors_inf(1,:), tb_vectors_inf(2,:), tb_vectors_inf(3,:), 'Color', 'b')
plot3([0 mtb_vectors_inf(1,:)], [0 mtb_vectors_inf(2,:)], [0 mtb_vectors_inf(3,:)], 'Color','b')
hold off

%% -----------------------------------------------------------------------
% OPTIMIZATION
% ------------------------------------------------------------------------
% Use the BP_solver function to optimize the normal vectors
[normal, exitflag] = BP_solver(tb_vectors_ap, tb_vectors_ml, [normal_ap0; normal_ml0;0;0]);
normal_ap = normal(1:3,1)/norm(normal(1:3,1));
normal_ml = normal(4:6,1)/norm(normal(4:6,1));

% Display results and dot product of the optimized normals
disp('      AP normal')
disp('   Initial    new')
disp([normal_ap0 normal_ap])
disp('      ML normal')
disp('   Initial    new')
disp([normal_ml0 normal_ml])
dot(normal_ap, normal_ml)

%% -----------------------------------------------------------------------
% VISUALIZATION OF BENDING PLANES
% ------------------------------------------------------------------------
% Define a point on the plane
point_on_plane = [0, 0, 0];

% Create a grid of points for the AP bending plane
[x_ap, y_ap] = meshgrid(linspace(-50, 30, 50), linspace(-10, 90, 50));
z_ap = -(normal_ap(1)*(x_ap - point_on_plane(1)) + normal_ap(2)*(y_ap - point_on_plane(2))) / normal_ap(3) + point_on_plane(3);

% Create a grid of points for the ML bending plane
[z_ml, y_ml] = meshgrid(linspace(-40, 30, 50), linspace(-10, 90, 50));
x_ml = -(normal_ml(3)*(z_ml - point_on_plane(3)) + normal_ml(2)*(y_ml - point_on_plane(2))) / normal_ml(1) + point_on_plane(1);

% Visualization of bending planes
figure(1)
view(-30,45)
camup([-1 0 0])
ax = gca;
ax.YRuler.FirstCrossoverValue  = 30; % X crossover with Y axis
ax.YRuler.SecondCrossoverValue  = 40;
ax.ZRuler.FirstCrossoverValue  = -20; % X crossover with Y axis
ax.ZRuler.SecondCrossoverValue  = 30;
hold on
f(3) = surf(x_ap, y_ap, z_ap, 'FaceAlpha', 0.2, 'EdgeColor','none', 'FaceColor','r','DisplayName', 'AP bending Plane');
f(4) = surf(x_ml, y_ml, z_ml, 'FaceAlpha', 0.2, 'EdgeColor','none','FaceColor','b','DisplayName', 'ML bending Plane');
legend(f)
box on
ax.BoxStyle = 'full';
xlabel('x')
ylabel('y')
zlabel('z')
hold off

% Save the optimized normals
save('normal_ap.mat', "normal_ap");
save('normal_ml.mat', "normal_ml");
