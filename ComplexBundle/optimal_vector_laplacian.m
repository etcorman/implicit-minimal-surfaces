function [W,Ae] = optimal_vector_laplacian(Src, conn, ifcplx, s)
% Integration formula from: "Globally Optimal Direction Fields" Section 6.1.1

if ~exist('ifcplx','var')
    ifcplx = false;
end
if ~exist('s','var')
    s = 0;
end

nv = Src.nv;
nf = Src.nf;
K = conn.K_tri;

idI = [1 2 3 1 2 3 2 3 1];
idJ = [1 2 3 2 3 1 1 2 3];
I = Src.T(:,idI);
J = Src.T(:,idJ);

ide  = abs(Src.T2E);
ides = sign(Src.T2E);
rho = exp(-1i*ides.*conn.para_trans_v2v(ide));

idK = abs(K) < 0.1;
Ki = K(idK);

if nargout > 1 || s ~= 0
    Oae = 2*(exp(1i*K) - 1 - 1i*K + K.^2/2 + 1i*K.^3/6)./(K.^4);
    Oae(idK) = (Ki.^5*1i)/181440 + Ki.^4/20160 - (Ki.^3*1i)/2520 - Ki.^2/360 + (Ki*1i)/60 + 1/12;
    Saeo = rho.*Oae;
    Sae = Src.area.*[ones(nf,3)/6, Saeo, conj(Saeo)];
    Ae = sparse(I, J, Sae, nv, nv);
    
    if ~ifcplx % real representation of vectors
        Ae = [real(Ae),-imag(Ae); imag(Ae), real(Ae)];
    end
    
    Ae = (Ae + Ae')/2;
end

f1K = (3 + 1i*K + K.^4/24 - 1i*K.^5/60 + exp(1i*K).*(-3 + 2*1i*K + K.^2/2))./K.^4;
f1K(idK) = - (Ki.^5*1i)/17280 - Ki.^4/2688 + (Ki.^3*1i)/504 + Ki.^2/120 - (Ki*1i)/24;
f2K = (4 + 1i*K - 1i*K.^3/6 - K.^4/12 + 1i*K.^5/30 + exp(1i*K).*(-4 + 3*1i*K + K.^2))./K.^4;
f2K(idK) =- (Ki.^5*7i)/51840 - Ki.^4/1120 + (Ki.^3*5i)/1008 + Ki.^2/45 - (Ki*1i)/24 - 1/4;

Swd = [(Src.cot_corner_angle(:,2)+Src.cot_corner_angle(:,3) + K.^2.*(3*Src.cot_corner_angle(:,1)+Src.cot_corner_angle(:,2)+Src.cot_corner_angle(:,3))/90), ...
       (Src.cot_corner_angle(:,3)+Src.cot_corner_angle(:,1) + K.^2.*(3*Src.cot_corner_angle(:,2)+Src.cot_corner_angle(:,3)+Src.cot_corner_angle(:,1))/90), ...
       (Src.cot_corner_angle(:,1)+Src.cot_corner_angle(:,2) + K.^2.*(3*Src.cot_corner_angle(:,3)+Src.cot_corner_angle(:,1)+Src.cot_corner_angle(:,2))/90)]/2;
Swo = [((Src.cot_corner_angle(:,1)+Src.cot_corner_angle(:,2)+2*Src.cot_corner_angle(:,3)).*f1K + Src.cot_corner_angle(:,3).*f2K), ... 
       ((Src.cot_corner_angle(:,2)+Src.cot_corner_angle(:,3)+2*Src.cot_corner_angle(:,1)).*f1K + Src.cot_corner_angle(:,1).*f2K), ...
       ((Src.cot_corner_angle(:,3)+Src.cot_corner_angle(:,1)+2*Src.cot_corner_angle(:,2)).*f1K + Src.cot_corner_angle(:,2).*f2K)]*2;
Swo = Swo.*rho;

if s == 0
    Sw = [Swd, Swo, conj(Swo)];
    W = sparse(I, J, Sw, nv, nv);
else
    Swo = Swo - s*(K.*Saeo - 1i*rho/2);
    Sw = [Swd - s*K/6, Swo, conj(Swo)];
    W = sparse(I, J, Sw, nv, nv);
end

% Vector Laplacian
if ~ifcplx % real representation of vectors
    W  = [real(W),-imag(W); imag(W), real(W)];
end

W = (W + W')/2;
