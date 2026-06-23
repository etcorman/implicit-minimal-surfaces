function [X,T,idvx_set_bnd] = close_hole(Src, if_double_cover)
if ~exist('if_double_cover','var')
    if_double_cover = false;
end

idvx_set_bnd = Src.idx_bound_cell;
X = Src.X;
T = Src.T;

% Remesh triangle at bnd
ide = all(ismember(Src.E2V, Src.idx_bound), 2) & all(Src.E2T(:,1:2) ~= 0, 2);
if any(ide)
    error('Need boundary remeshing.');
end

if if_double_cover
    nc = length(Src.idx_bound);
    nv = 2*Src.nv - nc;
    ind = setdiff(1:2*Src.nv, Src.nv+Src.idx_bound);
    ind_inv = zeros(2*Src.nv,1);
    ind_inv(Src.nv+Src.idx_bound) = Src.idx_bound;
    ind_inv(ind) = 1:nv;

    X = [X; X];
    X = X(ind,:);
    T = ind_inv([T; Src.nv + T(:,[2 1 3])]);
else
    nv = Src.nv;
    nc = length(idvx_set_bnd);
    for i = 1:nc
        Xm = mean(Src.X(idvx_set_bnd{i},:));
        Tm = [idvx_set_bnd{i}(1:end-1), (nv+1)*ones(length(idvx_set_bnd{i})-1,1), idvx_set_bnd{i}(2:end); ...
            idvx_set_bnd{i}(end), nv+1, idvx_set_bnd{i}(1)];

        X = [X; Xm];
        T = [T; Tm];
        nv = nv + 1;
    end
end