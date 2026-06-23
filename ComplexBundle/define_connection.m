function connection = define_connection(M, dec, idtri_sing)
% Discrete Connection and Covariant Derivative for Vector Field Analysis and Design

chi = M.nf - M.ne + M.nv;

if ~exist('idtri_sing','var')
    idtri_sing = 1;
end

% Compute parallel transport from vertex to vertex
theta = M.corner_angle;
ang_sum = accumarray(M.T(:), theta(:));
K_vx = 2*pi - ang_sum;
K_vx(M.idx_bound) = K_vx(M.idx_bound) - pi;

% Rescale angle
rescale = 2*pi./ang_sum;
rescale(M.idx_bound) = rescale(M.idx_bound)/2;
new_angle = theta.*rescale(M.T);

% Triangle Gaussian curvature
K_tri = sum(new_angle,2) - pi;
assert(norm(sum(K_tri) - sum(K_vx)) < 1e-6, 'Gaussian curvature incompatible with angle defect.');

% Connection complex line bundle
if abs(chi) ~= 0
    % Unique singularity
    sing_ghost = zeros(M.nf,1);
    sing_ghost(idtri_sing) = chi;

    % Boundary and cycle constraints
    Ibound = sparse(1:length(M.edge_bound), M.edge_bound, ones(length(M.edge_bound),1), length(M.edge_bound), M.ne);
    bbound = zeros(length(M.edge_bound),1);
    A = [dec.d1p; Ibound];
    b =-[K_tri - 2*pi*sing_ghost; bbound];

    % Solve system
    om = quadprog(speye(M.ne,M.ne), zeros(M.ne,1), [],[], A, b);

    % Check constraint
    sing = (K_tri + dec.d1p*om)/(2*pi);
    assert(max(abs(sing - sing_ghost) < 1e-6), 'Singularity placement failure.');

    % Connection wrt new basis
    para_trans_v2v =-om;
    assert(max(abs(wrapToPi(dec.d1p*para_trans_v2v - K_tri))) < 1e-6, 'Gaussian curvature incompatible with angle defect.');
    if max(abs(wrapToPi(dec.d1p*para_trans_v2v)/chi - wrapToPi(K_tri)/chi)) > 1e-6
        warning('The new basis does not allow dividing by chi.');
    end
else
    error('Connection undefined for tori.')
end

% Structure connection
p = 1/abs(chi);
connection.para_trans_v2v = p*para_trans_v2v;
connection.K_tri = p*K_tri;
