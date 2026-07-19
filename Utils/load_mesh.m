function [M,dec,conn,const] = load_mesh(name)

[X,T] = readOBJ([name, '.obj']);
M = MeshInfo(X, T, true);

disp('Load landmarks...');
[idvx,idvx_set] = read_landmarks([name, '.pinned']);
disp(['  ', num2str(length(idvx)), ' points -- ', num2str(length(idvx_set)), ' curves.']);
if ~isempty(M.idx_bound)
    disp('Filling holes...');
    [X,T,idvx_set_bnd] = close_hole(M);
    M = MeshInfo(X, T, false);
    idvx_set = [idvx_set; idvx_set_bnd];
    disp(' holes filled.');
end
disp('Computing connection...');
dec = dec_tri(M);
conn = define_connection(M, dec);

const.idvx = idvx;
const.idvx_set = idvx_set;

% Check input correct
check_input(M);