function C = getSymplecticBasis(M)
%GETSYMPLECTICBASIS Attempts to find the symplectic eigenbasis C generated
%from the eigensystem of a matrix (such as a monodromy matrix) M.

[eigvecs,eigvals]=eig(M,'vector');

%We sort eigvals based on real part; we then
%use this sorting arrangement to sort corresponding
%eigenvectors. The reason that we sort is to ensure that the eigenvalues
%are always in the order [lambda^-1 a-bi a+bi lambda] where lambda is the
%real eigenvalue > 1 and a and b are the components of the complex
%eigenvalue pair.
[eigvals,sortorder] = sort(eigvals);
eigvecs = eigvecs(:,sortorder);

%To ensure consistency, we also multiply each eigenvector by the sign of
%its first entry so that the sign of the first entry will always be
%positive. I'm not sure this step is actually necessary, but it is
%reassuring.

eigvecs = eigvecs .* sign(eigvecs(1,:))

%We also exchange the two complex eigenvectors. Sometimes, the less-than
%sign needs to be replaced with a greater-than sign or vice-versa.
if eigvecs(2,3) < 0
    eigvecs = eigvecs(:,[1 3 2 4]);
end

%We extract each generalized eigenvector. I'll be honest: Some choices here
%are somewhat based on experimentation.
usigma = -eigvecs(:,4);
usigmainv = eigvecs(:,1);
uomegap = real(eigvecs(:,3));
vomegap = imag(eigvecs(:,3));

%We order the eigenvectors into a "test" C.
C = [usigma uomegap usigmainv vomegap];

%Because C.' * J * C likely won't equal J yet, we find
%scalings for the columns of C by computing this expression.
Jproduct = C.'*Jmatrix(4)*C;

s1 = sqrt(Jproduct(1,3));
s2 = sqrt(Jproduct(2,4));

C = [usigma/s1 uomegap/s2 usigmainv/s1 vomegap/s2];

end

