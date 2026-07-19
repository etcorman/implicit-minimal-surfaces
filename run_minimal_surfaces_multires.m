% Implicit Minimal Surfaces for Bijective Correspondences
% Upscaling
clear all; close all; clc;
addpath(genpath(pwd));

%% Prepare file list
path_mesh = 'Mesh/';
name = 'hands_easy'; % sitting_armadillo
save_path = ['Results/', name, '/'];

if ~exist(save_path, 'dir')
    mkdir(save_path);
end

mesh_pathA = @(it) [path_mesh, name, '/A_' num2str(it), '.obj'];
mesh_pathB = @(it) [path_mesh, name, '/B_' num2str(it), '.obj'];

% Load landmarks
[X1_ref,T1_ref] = readOBJ([path_mesh, name, '/A.obj']);
M1_ref = MeshInfo(X1_ref, T1_ref, true);
[idvx1_ref,idvx_set1_ref] = read_landmarks([path_mesh, name, '/A.pinned']);

[X2_ref,T2_ref] = readOBJ([path_mesh, name, '/B.obj']);
M2_ref = MeshInfo(X2_ref, T2_ref, true);
[idvx2_ref,idvx_set2_ref] = read_landmarks([path_mesh, name, '/B.pinned']);

%% Main loop
it_max = 2;

for it = 1:it_max

    %% Upscale M1 and M2
    % If not the first iteration save old info
    if it > 1
        M1i = M1;
        conn1i = conn1;
        tri_sing1i = tri_sing1;

        M2i = M2;
        conn2i = conn2;
        tri_sing2i = tri_sing2;

        z_old = z;
    end

    % Down sample meshes
    if it < it_max
        [X1,T1] = readOBJ(mesh_pathA(it));
        [X2,T2] = readOBJ(mesh_pathB(it));
    else
        [X1,T1] = readOBJ([path_mesh, name, '/A.obj']);
        [X2,T2] = readOBJ([path_mesh, name, '/B.obj']);
    end

    % Load mesh 1
    M1 = MeshInfo(X1, T1, true);
    dec1 = dec_tri(M1);

    if it == 1
        tri_sing1 = 1;
    else
        m = (M1i.X(M1i.T(tri_sing1i,1),:) + M1i.X(M1i.T(tri_sing1i,2),:) + M1i.X(M1i.T(tri_sing1i,3),:))/3;
        bar  = (M1.X(M1.T(:,1),:) + M1.X(M1.T(:,2),:) + M1.X(M1.T(:,3),:))/3;
        [~,tri_sing1] = min(sum((bar - m).^2,2));
    end
    conn1 = define_connection(M1, dec1, tri_sing1);

    % Load mesh 2
    M2 = MeshInfo(X2, T2, true);
    dec2 = dec_tri(M2);

    if it == 1
        tri_sing2 = 1;
    else
        m = (M2i.X(M2i.T(tri_sing2i,1),:) + M2i.X(M2i.T(tri_sing2i,2),:) + M2i.X(M2i.T(tri_sing2i,3),:))/3;
        bar  = (M2.X(M2.T(:,1),:) + M2.X(M2.T(:,2),:) + M2.X(M2.T(:,3),:))/3;
        [~,tri_sing2] = min(sum((bar - m).^2,2));
    end
    conn2 = define_connection(M2, dec2, tri_sing2);

    % Projection directions
    if it == 1
        [~,id] = sort(max(M1.X) - min(M1.X));
        e1 = zeros(1,3); e1(id(3)) = 1;
        e2 = zeros(1,3); e2(id(2)) = 1;
    end

    %% Define Laplacian
    [W1,Ae1] = optimal_vector_laplacian(M1, conn1, true);
    [z1,l1] = eigs(W1, Ae1, 1, 'sm');
    
    [W2,Ae2] = optimal_vector_laplacian(M2, conn2, true);
    [z2,l2] = eigs(W2, Ae2, 1, 'sm');
    
    l = l1 + l2;

    %% Pinning function
    % Distance to constraints
    idvx1 = unique(knnsearch(M1.X, M1_ref.X(idvx1_ref,:)), 'stable');
    idvx2 = unique(knnsearch(M2.X, M2_ref.X(idvx2_ref,:)), 'stable');

    % Show distances
    D2 = [];
    for i = 1:length(idvx1)
        D2i = square_distances_4D(idvx1(i), idvx2(i), M1, dec1, M2, dec2, false);
        if isempty(D2)
            D2 = D2i;
        else
            D2 = min(D2, D2i);
        end
    end

    sigma_pin = 2;
    f_pin = 1 - exp(-D2*sigma_pin);
    if isempty(f_pin)
        f_pin = 1;
    end

    %% Compute init vector field
    if it > 1
        zi = transfer_complex_field(z_old, M1i, conn1i, M1, M2i, conn2i, M2);
    else
        % Nearest point map
        v2f = mesh_nearest_point(M1.X, M2);
        f2v = mesh_nearest_point(M2.X, M1);

        % Smooth vector field
        zi = vertex2face_initialization(v2f, f2v, M1, conn1, dec1, M2, conn2, dec2, 3000);
    end

    % Save init
    save_v2v_transfer(save_path, ['_init_it', num2str(it)], zi, M1, conn1, M2, conn2, e1, e2);

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
    writeObj([save_path, 'A_inter_it', num2str(it), '.obj'], X1_tot, T_tot, UV1, T_tot);
    writeObj([save_path, 'B_inter_it', num2str(it), '.obj'], X2_tot, T_tot, UV1, T_tot);

    T_tri_col = repmat((1:size(T_tot,1))', [1,3]);
    writeObj([save_path, 'A_inter_tri_it', num2str(it), '.obj'], X1_tot, T_tot, tri_col, T_tri_col, [], [], E2V_2(:,2:3));
    writeObj([save_path, 'B_inter_tri_it', num2str(it), '.obj'], X2_tot, T_tot, tri_col, T_tri_col, [], [], E2V_1(:,2:3));

end
