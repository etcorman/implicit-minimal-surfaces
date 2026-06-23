function [T_tot,X1_tot,X2_tot,tri_col,E2V_1,E2V_2,nnz] = compute_intersection_mesh(z, M1, conn1, dec1, M2, conn2, dec2, if_safe_mode)

if ~exist('if_safe_mode', 'var')
    if_safe_mode = true;
end
tol = 1e-8;

[omega_EV,omega_VE,sing_EE,sing_FV,sing_VF] = singularity_matrices(z, dec1, conn1, dec2, conn2);

% Compute edge-edge intersection
[X_EE1,X_EE2,Inter_index,E2V_1,E2V_2] = locate_edge_edge_singularities(z, M1, M2, omega_EV, omega_VE, sing_EE, tol);

% Vertex positons
[X2_1,npt1] = locate_vertex_face_singularities(M1, conn1, z, tol);
[X1_2,npt2] = locate_vertex_face_singularities(M2, conn2, z.', tol);

X1_tot = [M1.X; X2_1; X_EE1];
X2_tot = [X1_2; M2.X; X_EE2];
nnz = sum(npt1 ~= 1) + sum(npt2 ~= 1);

% Triangles of intersection mesh
T_tot = zeros(5*M1.nf+5*M2.nf,3);
tri_col = zeros(size(T_tot,1),2);
nf = 0;
for i = 1:M1.nf
    ide1 = abs(M1.T2E(i,:));
    for j = 1:M2.nf
        ide2 = abs(M2.T2E(j,:));

        inter_edge = Inter_index(ide1,ide2);
        id1 = sing_VF(M1.T(i,:),j) ~= 0;
        id2 = sing_FV(i,M2.T(j,:)) ~= 0;
        if any(id1) || any(id2) || any(inter_edge ~= 0, 'all')
            % List of vertices in the polygon
            idvx_list = [M1.T(i,id1)'; M1.nv + M2.T(j,id2)'; double(setdiff(inter_edge, 0))];

            % Sort vertices by angle on M1
            Xij1 = X1_tot(idvx_list,:) - mean(X1_tot(idvx_list,:));
            Xij1 = Xij1./sqrt(sum(Xij1.^2,2));
            if if_safe_mode
                assert(all(sum(Xij1.*M1.normal(i,:),2) < tol), 'Intersection not tangent to triangle.');
            end
            angij = atan2(sum(Xij1.*M1.e2r(i,:),2), sum(Xij1.*M1.e1r(i,:),2));
            [~,id] = sort(angij);
            idvx_list = idvx_list(id);
            
            % Triangulate polygons
            Tri_polygon = [idvx_list(1)*ones(length(idvx_list)-2,1), idvx_list(2:end-1), idvx_list(3:end)];
            id = nf+1:nf+size(Tri_polygon,1);
            T_tot(id,:) = Tri_polygon;
            tri_col(id,1) = i;
            tri_col(id,2) = j;
            nf = nf + size(Tri_polygon,1);
        end
    end
end
T_tot = T_tot(1:nf,:);
tri_col = tri_col(1:nf,:);
end

