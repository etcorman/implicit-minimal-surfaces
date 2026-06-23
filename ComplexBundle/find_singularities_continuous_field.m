function [pt,idtri_sing,l] = find_singularities_continuous_field(Src, conn, z, tol)
% Compute singularities using nonlinear interpolation in triangles:
% "Discrete Connection and Covariant Derivative for Vector Field Analysis and Design"

if ~exist('tol','var')
    tol = 1e-8;
end

% Compute singularities
omega = conn.para_trans_v2v;
rot = sign(Src.T2E).*omega(abs(Src.T2E));
K = conn.K_tri;

% Compute rotation
z_tri = z(Src.T);
ang_tri = angle(z_tri);
om = wrapToPi(ang_tri(:,[2 3 1]) - ang_tri(:,[1 2 3]) - rot);
tot_rot = sum(om,2) + K;
idtri_sing = find(abs(tot_rot) > pi/4);
index = round(tot_rot(idtri_sing)/(2*pi));
if ~isempty(idtri_sing)
    assert(max(abs(2*pi*index - tot_rot(idtri_sing))) < 1e-6, 'Non-quantized singu.');
else
    warning('No singularity found...');
end

% Compute singularities
pt = zeros(Src.nf,3);
l = ones(Src.nf,3)/3;
for tri = idtri_sing'
    zn = abs(z_tri(tri,:));
    zn = zn/max(zn);
    Ki = K(tri);

    a = Ki/3 + (-2*om(tri,:) + om(tri,[2 3 1]) + om(tri,[3 1 2]))/3;
    le = zn(:,[2 3 1]).*zn(:,[3 1 2]).*sin(om(tri,[2 3 1]) + a(:,[2 3 1]));
    le = le/sum(le);
    for t = linspace(0,1,10)
        znt = zn + 1 - t;
        for j = 1:20
            eqsin = sin(om(tri,[2 3 1]) + le*Ki*t + (1-t)*a(:,[2 3 1]));
            eqcos = cos(om(tri,[2 3 1]) + le*Ki*t + (1-t)*a(:,[2 3 1]));
            A = [znt(1), znt(2)*eqcos(3), znt(3)*eqcos(2); 0, znt(2)*eqsin(3),-znt(3)*eqsin(2); 1, 1, 1];
            F = A*le' - [0;0;1];
            Jf = A + Ki*t*[0,-znt(3)*eqsin(2)*le(3),-znt(2)*eqsin(3)*le(2); 0,-znt(3)*eqcos(2)*le(3), znt(2)*eqcos(3)*le(2); 0, 0, 0];

            le = le - (Jf\F)';

            if max(abs(F)) < tol
                break;
            end
        end
    end
    l(tri,:) = le;
    pt(tri,:) = l(tri,1)*Src.X(Src.T(tri,1),:) + l(tri,2)*Src.X(Src.T(tri,2),:) + l(tri,3)*Src.X(Src.T(tri,3),:);

    % Check if 0 found
    if  max(abs(F)) >= tol
        warning('Zero not found.');
    end
end
pt = pt(idtri_sing,:);

if nargout >= 3
    l = l(idtri_sing,:);
end