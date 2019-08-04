function  resultObj = ctranspose(D)

D.adjoint = xor(D.adjoint,1);
resultObj = D;
