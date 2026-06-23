function save_v2v_transfer(save_path, suffix, z, M1, conn1, M2, conn2, e1, e2, idvx1, idvx2)

tol = 1e-10;

% Param
UV1 = [sum(M1.X.*e1,2), sum(M1.X.*e2,2)];
UV1 = (UV1 - min(UV1))/max((max(UV1) - min(UV1)));
writeObj([save_path, 'A', suffix, '.obj'], M1.X, M1.T, UV1, M1.T);

% Compute correspondences
[~,idtri,l_tri] = locate_vertex_face_singularities(M1, conn1, z, tol);
UV2 = UV1(M1.T(idtri,1),:).*l_tri(:,1) + UV1(M1.T(idtri,2),:).*l_tri(:,2) + UV1(M1.T(idtri,3),:).*l_tri(:,3);

% Save OBJs
writeObj([save_path, 'B', suffix, '.obj'], M2.X, M2.T, UV2, M2.T);

% Export constraints
if exist('idvx1','var') && ~isempty(idvx1)
    writeObj([save_path, 'A_lm.obj'], M1.X(idvx1,:), []);
    writeObj([save_path, 'B_lm.obj'], M2.X(idvx2,:), []);
end
end

