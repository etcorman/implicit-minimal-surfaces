function [omega_EV,omega_VE,sing_EE,sing_FV,sing_VF] = singularity_matrices(z, dec1, conn1, dec2, conn2)

tol = 1e-8;

% Field rotations
omega_EV = wrapToPi(dec1.d0p*angle(z) - conn1.para_trans_v2v);
omega_VE = wrapToPi(angle(z)*dec2.d0p' - conn2.para_trans_v2v');

% Singularity matrices
sing_EE = (dec1.d0p*omega_VE - omega_EV*dec2.d0p')/(2*pi);
sing_FV = (dec1.d1p*omega_EV + conn1.K_tri)/(2*pi);
sing_VF = (omega_VE*dec2.d1p' + conn2.K_tri')/(2*pi);

assert(max(abs(vec(sing_FV - round(sing_FV)))) < tol);
assert(max(abs(vec(sing_VF - round(sing_VF)))) < tol);
assert(max(abs(vec(sing_EE - round(sing_EE)))) < tol);

sing_EE = cast(round(sing_EE), 'int8');
sing_FV = cast(round(sing_FV), 'int8');
sing_VF = cast(round(sing_VF), 'int8');