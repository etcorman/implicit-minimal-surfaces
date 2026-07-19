% Implicit Minimal Surfaces for Bijective Correspondences
clear all; close all; clc;
addpath(genpath(pwd));

%% Prepare file list
path_mesh = 'Mesh/';
name = 'hands_easy'; % hands_easy teddy sitting_armadillo
save_path = ['Results/', name, '/'];

if ~exist(save_path, 'dir')
    mkdir(save_path);
end

%% Load meshes
disp('Load mesh 1...');
[M1,dec1,conn1,const1] = load_mesh([path_mesh, name, '/A']);
disp('done loading mesh 1.');

disp('Load mesh 2...');
[M2,dec2,conn2,const2] = load_mesh([path_mesh, name, '/B']);
disp('done loading mesh 2.');

% Best projection
[~,id] = sort(max(M1.X) - min(M1.X));
e1 = zeros(1,3); e1(id(3)) = 1;
e2 = zeros(1,3); e2(id(2)) = 1;

%% Define Laplacian
[W1,Ae1] = optimal_vector_laplacian(M1, conn1, true);
[z1,l1] = eigs(W1, Ae1, 1, 'sm');

[W2,Ae2] = optimal_vector_laplacian(M2, conn2, true);
[z2,l2] = eigs(W2, Ae2, 1, 'sm');

l = l1 + l2;

%% Pinning function
% Distance to constraints
D2 = [];
for i = 1:length(const1.idvx)
    D2i = square_distances_4D(const1.idvx(i), const2.idvx(i), M1, dec1, M2, dec2, false);
    if isempty(D2)
        D2 = D2i;
    else
        D2 = min(D2, D2i);
    end
end

for i = 1:length(const1.idvx_set)
    D2i = square_distances_4D(const1.idvx_set{i}, const2.idvx_set{i}, M1, dec1, M2, dec2, false);
    if isempty(D2)
        D2 = D2i;
    else
        D2 = min(D2, D2i);
    end
end

% Pinning function
sigma_pin = 2;
f_pin = 1 - exp(-D2*sigma_pin);
if isempty(f_pin)
    f_pin = 1;
end

%% Compute init vector field
% Nearest point map
v2f = mesh_nearest_point(M1.X, M2);
f2v = mesh_nearest_point(M2.X, M1);

% Smooth vector field
zi = vertex2face_initialization(v2f, f2v, M1, conn1, dec1, M2, conn2, dec2, 3000);

% Save init
save_v2v_transfer(save_path, '_init', zi, M1, conn1, M2, conn2, e1, e2);

%% LBFGS
opts = optimoptions("fminunc",Display="none", MaxFunctionEvaluations=Inf,MaxIterations=3000);
opts.Display = 'iter';
opts.HessianApproximation = 'lbfgs';
opts.SpecifyObjectiveGradient = true;
opts.OptimalityTolerance = 1e-5;

% Set parameters
xi = [real(zi(:)); imag(zi(:))];
lambda = 100*l;
fun = @(x) GinzburgLandau_lbfgs(x, lambda, W1, Ae1, dec1, W2, Ae2, dec2, f_pin);

% Actual optimization
[x,fval,exitflag,output] = fminunc(fun, xi, opts);

u = reshape(x(1:M1.nv*M2.nv), [M1.nv,M2.nv]);
v = reshape(x(M1.nv*M2.nv+1:2*M1.nv*M2.nv), [M1.nv,M2.nv]);
z = complex(u, v);

%% Plot vertex mapping
show_v2v_mapping(z, M1, conn1, M2, conn2);

%% Mesh intersection
[T_tot,X1_tot,X2_tot,tri_col,E2V_1,E2V_2] = compute_intersection_mesh(z, M1, conn1, dec1, M2, conn2, dec2);

% Export obj file
UV1 = [sum(X1_tot.*e1,2), sum(X1_tot.*e2,2)];
UV1 = (UV1 - min(UV1))/max((max(UV1) - min(UV1)));
writeObj([save_path, 'A_inter.obj'], X1_tot, T_tot, UV1, T_tot);
writeObj([save_path, 'B_inter.obj'], X2_tot, T_tot, UV1, T_tot);

T_tri_col = repmat((1:size(T_tot,1))', [1,3]);
writeObj([save_path, 'A_tri.obj'], X1_tot, T_tot, tri_col, T_tri_col, [], [], E2V_2(:,2:3));
writeObj([save_path, 'B_tri.obj'], X2_tot, T_tot, tri_col, T_tri_col, [], [], E2V_1(:,2:3));
