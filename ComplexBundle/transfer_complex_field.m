function z = transfer_complex_field(z_init, M1i, conn1i, M1, M2i, conn2i, M2)

% Project vertices of M1 to faces of M1i
[tri1,bar1] = mesh_nearest_point(M1.X, M1i);
[tri2,bar2] = mesh_nearest_point(M2.X, M2i);

% Interpolate complex field
Phi1 = interpolation_basis_eval(M1i, conn1i, tri1, bar1);
Phi2 = interpolation_basis_eval(M2i, conn2i, tri2, bar2);

z = zeros(M1.nv,M2.nv);
for i = 1:M1.nv
    zi = z_init(M1i.T(tri1(i),:),:);
    for j = 1:M2.nv
        zij = zi(:,M2i.T(tri2(j),:));

        z(i,j) = sum(Phi1(i,:).'.*zij.*Phi2(j,:), 'all');
    end
end
end

function Phi = interpolation_basis_eval(M, conn, tri, bar)
% Half-edge parallel transport
rho_half_egde = sign(M.T2E).*conn.para_trans_v2v(abs(M.T2E));
rot_vx2f = zeros(M.nf,3);
rot_vx2f(:,2) = rot_vx2f(:,1) + conn.K_tri/3 - rho_half_egde(:,1);
rot_vx2f(:,3) = rot_vx2f(:,2) + conn.K_tri/3 - rho_half_egde(:,2);
% err = rot_vx2f(:,3) + conn.K_tri/3 - rot_vx2f(:,1) - rho_half_egde(:,3);
% assert(max(abs(wrapToPi(err))) < 1e-6, 'Mismatched vertex to corner parallel transport');

ang =   rot_vx2f(tri,:) ...
      + bar(:,[2 3 1]).*conn.K_tri(tri)/3 ...
      - bar(:,[3 1 2]).*conn.K_tri(tri)/3;

Phi = bar.*exp(1i*ang);
end