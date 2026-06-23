function D2 = square_distances_4D(idvx1, idvx2, M1, dec1, M2, dec2, if_plot)

if ~exist('if_plot','var') || isempty(if_plot)
    if_plot = flase;
end

% Check input data
if size(idvx1,1) == 1 && size(idvx1,2) ~= 1
    idvx1 = idvx1';
end
if size(idvx2,1) == 1 && size(idvx2,2) ~= 1
    idvx2 = idvx2';
end
assert(size(idvx1,2) == 1 && size(idvx2,2) == 1, 'Not the same number of correspondences.');
assert(all(idvx1 > 0) && all(idvx1 <= M1.nv), 'Variable idvx1 is not a valid list of indices.');
assert(all(idvx2 > 0) && all(idvx2 <= M2.nv), 'Variable idvx2 is not a valid list of indices.');
assert(length(unique(idvx1)) == length(idvx1), 'Variable idvx1 define correspondences to the same point.');
assert(length(unique(idvx2)) == length(idvx2), 'Variable idvx2 define correspondences to the same point.');

% Geodesic in heat on M1
d1 = [];
for i = 1:length(idvx1)
    b = zeros(M1.nv,1);
    b(idvx1(i)) = 1;
    di = geodesic_heat(b, M1, dec1);
    if isempty(d1)
        d1 = di;
    else
        d1 = min(d1, di);
    end
end

% Geodesic in heat on M2
d2 = [];
for i = 1:length(idvx2)
    b = zeros(M2.nv,1);
    b(idvx2(i)) = 1;
    di = geodesic_heat(b, M2, dec2);
    if isempty(d2)
        d2 = di;
    else
        d2 = min(d2, di);
    end
end

% Squared distance in 4D
D2 = d1.^2 + d2.^2';

% Plot stuff
if if_plot
    figure(1);
    subplot(1,2,1);
    trisurf(M1.T, M1.X(:,1), M1.X(:,2), M1.X(:,3), d1, 'facecolor', 'interp');
    hold on; scatter3(M1.X(idvx1,1), M1.X(idvx1,2), M1.X(idvx1,3), 100, 'filled'); hold off;
    axis equal; colorbar;
    subplot(1,2,2);
    trisurf(M2.T, M2.X(:,1), M2.X(:,2), M2.X(:,3), d2, 'facecolor', 'interp');
    hold on; scatter3(M2.X(idvx2,1), M2.X(idvx2,2), M2.X(idvx2,3), 100, 'filled'); hold off;
    axis equal; colorbar;
    pause
end
end

function d = geodesic_heat(rhs, Src, dec)
% Heat diffusion
dt = mean(sqrt(Src.SqEdgeLength))^2;
hd = (dt*dec.W + dec.star0p)\rhs;

% Renormalization
[gradf1,gradf2,gradf3] = gradient_fem(hd , Src.X, Src.T, Src.area);
norm_grad = sqrt(gradf1.^2 + gradf2.^2 + gradf3.^2);
gradf1 = gradf1./norm_grad; gradf2 = gradf2./norm_grad; gradf3 = gradf3./norm_grad;

% Divergence
div_vf = divergence_fem(gradf1, gradf2, gradf3, Src.X, Src.T);

% Compute distance
d =-(dec.W + 1e-7*dec.star0p)\div_vf;
d = d - min(d(:));
end

function [df1,df2,df3] = gradient_fem(f, V, T, A)
nf = size(T, 1);
n = size(f,2);

if ~exist('A', 'var')
    normal = cross(V(T(:,1),:) - V(T(:,2),:), V(T(:,1),:) - V(T(:,3),:));
    A = sqrt(sum(normal.^2, 2))/2;
end
N = cross( V(T(:,1),:) - V(T(:,2),:), V(T(:,1),:) - V(T(:,3),:));
N = N./repmat(sqrt(sum(N.^2, 2)), [1, 3]);

idj = [2 3 1];
idK = [3 1 2];
df1 = zeros(nf, n);
df2 = zeros(nf, n);
df3 = zeros(nf, n);
for i = 1:3
    j = idj(i);
    k = idK(i);
    eij = cross(N, V(T(:,j),:) - V(T(:,i),:), 2)./(2*A);
    df1 = df1 + f(T(:,k),:) .* eij(:,1);
    df2 = df2 + f(T(:,k),:) .* eij(:,2);
    df3 = df3 + f(T(:,k),:) .* eij(:,3);
end

end

function [div_vf] = divergence_fem(vf1, vf2, vf3, V, T)
nf = size(T, 1);
nv = size(V, 1);
n = size(vf1,2);

N = cross( V(T(:,1),:) - V(T(:,2),:), V(T(:,1),:) - V(T(:,3),:));
N = N./repmat(sqrt(sum(N.^2, 2)), [1, 3]);

idj = [2 3 1];
idK = [3 1 2];
div_vf = zeros(nv, n);
for i = 1:3
    j = idj(i);
    k = idK(i);
    eij = cross(N, V(T(:,j),:) - V(T(:,i),:), 2);
    scalar_prod = vf1 .* eij(:,1) + vf2 .* eij(:,2) + vf3 .* eij(:,3);

    Tk = sparse((1:nf)', T(:,k), 1, nf, nv);
    div_vf = div_vf + Tk'*scalar_prod;
end
div_vf = 0.5 * div_vf;
end