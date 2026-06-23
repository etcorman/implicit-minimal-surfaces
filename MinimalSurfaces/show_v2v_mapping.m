function show_v2v_mapping(z, M1, conn1, M2, conn2)

X2_1 = locate_vertex_face_singularities(M1, conn1, z);
X1_2 = locate_vertex_face_singularities(M2, conn2, z.');

% Draw
col_1 = M1.X; col_1 = (col_1 - min(col_1))./(max(col_1) - min(col_1));
col_2 = M2.X; col_2 = (col_2 - min(col_2))./(max(col_2) - min(col_2));
figure;
subplot(2,2,1);
trisurf(M1.T, M1.X(:,1), M1.X(:,2), M1.X(:,3), 1);
hold on; scatter3(M1.X(:,1), M1.X(:,2), M1.X(:,3), 100, col_1, 'filled'); hold off;
axis equal;
title('A');
subplot(2,2,2);
trisurf(M1.T, X1_2(:,1), X1_2(:,2), X1_2(:,3), 1);
hold on; scatter3(X1_2(:,1), X1_2(:,2), X1_2(:,3), 100, col_1, 'filled'); hold off;
axis equal;
title('A to B');

subplot(2,2,3);
trisurf(M2.T, M2.X(:,1), M2.X(:,2), M2.X(:,3), 1);
hold on; scatter3(M2.X(:,1), M2.X(:,2), M2.X(:,3), 100, col_2, 'filled'); hold off;
axis equal;
title('B');
subplot(2,2,4);
trisurf(M2.T, X2_1(:,1), X2_1(:,2), X2_1(:,3), 1);
hold on; scatter3(X2_1(:,1), X2_1(:,2), X2_1(:,3), 100, col_2, 'filled'); hold off;
title('B to A');
axis equal;
