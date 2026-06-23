function [f,df] = GinzburgLandau_lbfgs(x, lambda, W1, A1, dec1, W2, A2, dec2, f_pin)
n = size(W1,1);
m = size(W2,1);
X = reshape(x, [n,2*m]);
u = X(:,1:m);
v = X(:,m+1:2*m);

if ~exist('f_pin','var') || isempty(f_pin)
    f_pin = 1;
end

W1r = real(W1); W1i = imag(W1);
A1r = real(A1); A1i = imag(A1);
W2r = real(W2); W2i = imag(W2);
A2r = real(A2); A2i = imag(A2);

A = (W1r*u - W1i*v);
B = (W1r*v + W1i*u);

C = (A1r*u - A1i*v);
D = (A1r*v + A1i*u);

df_11 = A*A2r' - B*A2i';
df_12 = A*A2i' + B*A2r';
df_21 = C*W2r' - D*W2i';
df_22 = C*W2i' + D*W2r';

double_well = u.^2+v.^2 - f_pin;
f = sum(vec( u.*df_11 + v.*df_12 ))/2 ...
  + sum(vec( u.*df_21 + v.*df_22 ))/2 ... 
  + lambda*sum(vec( dec1.star0p*double_well.^2*dec2.star0p.' ))/4;
if nargout > 1
    double_well_weighted = lambda*double_well;
    df = vec([df_11 + df_21 + dec1.star0p*(double_well_weighted.*u)*dec2.star0p.', ...
              df_12 + df_22 + dec1.star0p*(double_well_weighted.*v)*dec2.star0p.']);
end
end