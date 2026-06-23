function [M,dec,conn] = load_mesh(name)

[X,T] = readOBJ([name, '.obj']);
M = MeshInfo(X, T, true);

disp('Load landmarks...');
[idvx1,idvx_set1] = read_landmarks([name, '.pinned']);
disp(['  ', num2str(length(idvx1)), ' points -- ', num2str(length(idvx_set1)), ' curves.']);
if ~isempty(M.idx_bound)
    disp('Filling holes...');
    [X,T,idvx_set_bnd1] = close_hole(M);
    M = MeshInfo(X, T, false);
    idvx_set1 = [idvx_set1; idvx_set_bnd1];
    disp(' holes filled.');
end
disp('Computing connection...');
dec = dec_tri(M);
conn = define_connection(M, dec);

check_input(M);