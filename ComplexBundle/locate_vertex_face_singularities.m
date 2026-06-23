function [X_src2tar,idtri,l_tri,npt,nout] = locate_vertex_face_singularities(Src, conn, Z, tol)
if ~exist('tol','var')
    tol = 1e-8;
end

idtri = zeros(size(Z,2),1);
l_tri = zeros(size(Z,2),3);
npt = zeros(size(Z,2),1);
nout = zeros(size(Z,2),1);
for j = 1:size(Z,2)
    [ptj,idtrj,l] = find_singularities_continuous_field(Src, conn, Z(:,j), tol);
    npt(j) = size(ptj,1);
    nout(j) = sum(any(l < 0,2) | any(l > 1,2));

    if length(idtrj) > 1
        normz = abs(Z(Src.T(idtrj,1),j)) + abs(Z(Src.T(idtrj,2),j)) + abs(Z(Src.T(idtrj,3),j));
        [~,k] = min(normz);
        idtrj = idtrj(k);
        l = l(k,:);
    elseif isempty(idtrj)
        idtri = 1;
        l = [1,0,0];
    end

    idtri(j) = idtrj;
    l_tri(j,:) = l;
end
X_src2tar = Src.X(Src.T(idtri,1),:).*l_tri(:,1) + Src.X(Src.T(idtri,2),:).*l_tri(:,2) + Src.X(Src.T(idtri,3),:).*l_tri(:,3);

disp(['Number of non-bijective field: ', num2str(sum(npt ~= 1))]);
disp(['Number of singu outside triangle: ', num2str(sum(nout))]);
end