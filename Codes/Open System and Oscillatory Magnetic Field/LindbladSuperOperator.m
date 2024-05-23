function L = LindbladSuperOperator(gamma,A,dim)

It = eye(dim);
L = gamma*(kron(conj(A),A)-0.5*kron(It,A'*A)-0.5*kron(A.'*conj(A),It));

end