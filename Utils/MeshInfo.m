function mesh = MeshInfo(X, T, if_rescale)

if ~exist('if_rescale', 'var')
	if_rescale = false;
end

if if_rescale
	area = sqrt(sum(cross(X(T(:,1),:) - X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:),2).^2,2))/2;
    area_tot = sum(area);
    bar = (X(T(:,1),:) + X(T(:,2),:) + X(T(:,3),:))/3;
    m = sum(bar.*area)/area_tot;
	X = (X - m)/sqrt(area_tot);
end

mesh.X = X;
mesh.T = T;
assert(size(T,2) == 3, 'Not a triangulations.');
[mesh.E2V, mesh.T2E, mesh.E2T, mesh.T2T] = connectivity(mesh.T);

mesh.nf = size(mesh.T,1);
mesh.nv = size(mesh.X,1);
mesh.ne = size(mesh.E2V,1);

% Normals and areas
mesh.normal = cross(mesh.X(mesh.T(:,1),:) - mesh.X(mesh.T(:,2),:), mesh.X(mesh.T(:,1),:) - mesh.X(mesh.T(:,3),:));
mesh.area = sqrt(sum(mesh.normal.^2, 2))/2;
mesh.normal = mesh.normal./repmat(sqrt(sum(mesh.normal.^2, 2)), [1, 3]);

A = sparse(mesh.T, repmat((1:mesh.nf)', [3,1]), repmat(mesh.area, [3,1]), mesh.nv, mesh.nf);
mesh.Nv = A*mesh.normal;
mesh.Nv = mesh.Nv./repmat(sqrt(sum(mesh.Nv.^2,2)), [1,3]);

% Edge length
mesh.SqEdgeLength = sum((mesh.X(mesh.E2V(:,1),:) - mesh.X(mesh.E2V(:,2),:)).^2, 2);

% Angles
mesh.corner_angle = angles_of_triangles(mesh.X, mesh.T);
mesh.cot_corner_angle = cot(mesh.corner_angle);

% Face local basis
mesh.e1r = mesh.X(mesh.T(:,1),:) - mesh.X(mesh.T(:,2),:);
mesh.e1r = mesh.e1r./sqrt(sum(mesh.e1r.^2,2));
mesh.e2r = cross(mesh.normal, mesh.e1r, 2);

% Extract boundaries
[idx_bound_cell,edge_bound_cell] = extract_all_boundary_curves(mesh.T, mesh.E2V);
mesh.idx_bound_cell = idx_bound_cell;
mesh.edge_bound_cell = edge_bound_cell;

tri_bound_cell = cell(size(edge_bound_cell));
for i = 1:length(tri_bound_cell)
    tri_bound_cell{i} = sum(mesh.E2T(edge_bound_cell{i},1:2), 2);
end
mesh.tri_bound_cell = tri_bound_cell;

mesh.idx_bound = cell2mat(idx_bound_cell);
mesh.edge_bound = cell2mat(edge_bound_cell);
mesh.tri_bound = cell2mat(tri_bound_cell);

mesh.idx_int = setdiff((1:mesh.nv)', mesh.idx_bound);
mesh.edge_int = setdiff((1:mesh.ne)', mesh.edge_bound);
mesh.tri_int = setdiff((1:mesh.nf)', mesh.tri_bound);
end
