function [z,zr,sing_FV,sing_VF] = vertex2face_initialization(v2f, f2v, M1, conn1, dec1, M2, conn2, dec2, it_cg)

if ~exist('it_cg','var') || isempty(it_cg)
    it_cg = 400;
end

sing = 2*pi;

eps = 1e-7;

om_EV = repmat(conn1.para_trans_v2v, [1,M2.nv]);
om_VE = repmat(conn2.para_trans_v2v', [M1.nv,1]);

assert(max(abs(vec( wrapToPi(dec1.d1p*om_EV - conn1.K_tri) ))) < 1e-6, 'Incompatible connection and Gaussian curvature.');
assert(max(abs(vec( wrapToPi(om_VE*dec2.d1p' - conn2.K_tri') ))) < 1e-6, 'Incompatible connection and Gaussian curvature.');

% Connection edge to vertex
FV_sing = sparse(f2v, 1:M2.nv, sing, M1.nf, M2.nv);
FV = FV_sing - conn1.K_tri;
assert(max(abs(sum(FV))) < eps);

W0d_1 = dec1.d0d'*dec1.star1d*dec1.d0d + eps*dec1.star0d; W0d_1 = (W0d_1 + W0d_1')/2;
G = W0d_1\FV;
omega_EV = dec1.star1d*dec1.d0d*G + om_EV;
assert(max(abs(vec(wrapToPi( dec1.d1p*omega_EV - FV_sing )))) < eps);

% Connection vertex to edge
VF_sing = sparse(1:M1.nv, v2f, sing, M1.nv, M2.nf);
VF = VF_sing - conn2.K_tri';
assert(max(abs(sum(VF, 2))) < eps);

W0d_2 = dec2.d0d'*dec2.star1d*dec2.d0d + eps*dec2.star0d; W0d_2 = (W0d_2 + W0d_2')/2;
G = W0d_2\VF';
omega_VE = (dec2.star1d*dec2.d0d*G)' + om_VE;
assert(max(abs(vec(wrapToPi( omega_VE*dec2.d1p' - VF_sing )))) < eps);

% Find smallest eigen vector
tensor_laplacian_mult = @(x) laplacian4d_mult(x, omega_EV, omega_VE, M1, dec1, M2, dec2);
% tensor_laplacian_mult = @(x) laplacian4d_mult2(x, omega_EV, FV_sing, omega_VE, VF_sing, M1, dec1, M2, dec2);
tensor_mass_mult = @(x) vec(dec1.star0p*reshape(x,[M1.nv,M2.nv])*dec2.star0p.');

% LOBPCG https://github.com/lobpcg/blopex/tree/master?tab=readme-ov-file
z = complex(randn(M1.nv*M2.nv,1), randn(M1.nv*M2.nv,1));
z = z/real(z'*tensor_mass_mult(z));
[z,mu] = lobpcg(z(:), tensor_laplacian_mult, tensor_mass_mult, [], 1e-3, it_cg, 1);
z = z*sqrt(M1.nv*M2.nv/real(z'*z));

z = reshape(z,[M1.nv,M2.nv]);
zr = [real(z), imag(z)];

sing_FV = (dec1.d1p*wrapToPi(dec1.d0p*angle(z)  - om_EV)  + conn1.K_tri  - FV_sing)/(2*pi);
sing_VF = (wrapToPi(angle(z)*dec2.d0p' - om_VE)*dec2.d1p' + conn2.K_tri' - VF_sing)/(2*pi);

end


function Wz = laplacian4d_mult(z, omega_EV, omega_VE, M1, dec1, M2, dec2)
z = reshape(z, [M1.nv,M2.nv]);
zEV = z(M1.E2V(:,1),:) - z(M1.E2V(:,2),:).*exp(1i*omega_EV);
zVE = z(:,M2.E2V(:,1)) - z(:,M2.E2V(:,2)).*exp(1i*omega_VE);
zEV_weighted = dec1.star1p*zEV*dec2.star0p';
zVE_weighted = dec1.star0p*zVE*dec2.star1p';

Int_E2V1_1 = sparse(M1.E2V(:,1), (1:M1.ne)', 1, M1.nv, M1.ne);
Int_E2V2_1 = sparse(M1.E2V(:,2), (1:M1.ne)', 1, M1.nv, M1.ne);
Int_E2V1_2 = sparse(M2.E2V(:,1), (1:M2.ne)', 1, M2.nv, M2.ne)';
Int_E2V2_2 = sparse(M2.E2V(:,2), (1:M2.ne)', 1, M2.nv, M2.ne)';

dzEV = Int_E2V1_1*zEV_weighted - Int_E2V2_1*(zEV_weighted.*exp(-1i*omega_EV));
dzVE = zVE_weighted*Int_E2V1_2 - (zVE_weighted.*exp(-1i*omega_VE))*Int_E2V2_2;

Wz = dzEV(:) + dzVE(:);
end

function Wz = laplacian4d_mult2(z, omega_EV, FV_sing, omega_VE, VF_sing, M1, dec1, M2, dec2)
z = reshape(z, [M1.nv,M2.nv]);

Wz1 = zeros(M1.nv,M2.nv);
for i = 1:M2.nv
    conn1.K_tri = FV_sing(:,i);
    conn1.para_trans_v2v = omega_EV(:,i);

    W = optimal_vector_laplacian_exact(M1, conn1, 1, true);
    Wz1(:,i) = W*z(:,i);
end
Wz1 = Wz1*dec2.star0p;

Wz2 = zeros(M1.nv,M2.nv);
for i = 1:M1.nv
    conn2.K_tri = VF_sing(i,:)';
    conn2.para_trans_v2v = omega_VE(i,:)';

    W = optimal_vector_laplacian_exact(M2, conn2, 1, true);
    Wz2(i,:) = z(i,:)*W.';
end
Wz2 = dec1.star0p*Wz2;

Wz = Wz1(:) + Wz2(:);
end