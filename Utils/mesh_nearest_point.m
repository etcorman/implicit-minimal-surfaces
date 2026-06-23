function [id_tri,l_tri] = mesh_nearest_point(Xs, varargin)
% Input: Xs : V x 3 query coordinates
%        X  : V x 3 target vertex coordinates
%        T  : V x 3 target triangle connectivity
% Input: Xs : V x 3 query coordinates
%        Src: structure containing X and T

if nargin == 2
    X = varargin{1}.X;
    T = varargin{1}.T;
elseif nargin == 3
    X = varargin{1};
    T = varargin{2};
else
    error('Wrong number of input arguments.')
end
nv = size(Xs,1);

% Resample target
X_sample = [X(T(:,1),:); ...
            X(T(:,2),:); ...
            X(T(:,3),:); ...
            (  X(T(:,1),:) +   X(T(:,2),:) + 0*X(T(:,3),:))/2; ...
            (  X(T(:,1),:) + 0*X(T(:,2),:) +   X(T(:,3),:))/2; ...
            (0*X(T(:,1),:) +   X(T(:,2),:) +   X(T(:,3),:))/2; ...
            (2*X(T(:,1),:) + 2*X(T(:,2),:) +   X(T(:,3),:))/5; ...
            (2*X(T(:,1),:) +   X(T(:,2),:) + 2*X(T(:,3),:))/5; ...
            (  X(T(:,1),:) + 2*X(T(:,2),:) + 2*X(T(:,3),:))/5; ...
            (2*X(T(:,1),:) +   X(T(:,2),:) +   X(T(:,3),:))/4; ...
            (  X(T(:,1),:) + 2*X(T(:,2),:) +   X(T(:,3),:))/4; ...
            (  X(T(:,1),:) +   X(T(:,2),:) + 2*X(T(:,3),:))/4];

id_tri = knnsearch(X_sample, Xs);
id_tri = mod(id_tri-1, size(T,1)) + 1;

% Compute nearest point on triangle
l_tri = zeros(nv,3);
e = ones(3,1);
for i = 1:nv
    y = Xs(i,:)';
    pt = X(T(id_tri(i),:),:)';
    A = [pt'*pt, e; e', 0];
    b = [pt'*y; 1];
    l = A\b;
    l_tri(i,:) = l(1:3);
    if any(isnan(l_tri(i,:)))
        l_tri(i,:) = [1, 0, 0];
    end
end
assert(max(abs( sum(l_tri,2) - 1 )) < 1e-6, 'Constraint failed');

% Reproject on nearest valid coordinates
id = any(l_tri <= 0,2);
l_tri(id,:) = max(l_tri(id,:), 0);
l_tri(id,:) = l_tri(id,:)./sum(l_tri(id,:), 2);
assert(all((l_tri(:) >= 0) & (l_tri(:) <= 1)));
