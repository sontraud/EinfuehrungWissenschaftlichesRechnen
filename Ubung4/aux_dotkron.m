function H = aux_dotkron(A,B)
%some sort of Kronecker tensor product
%used to assemble F'(c) fast in the inverse EIT algorithm

[m,n]=size(A);
[k,l]=size(B);

if n~=l
    error('dimension fault in dotkron');
end

%A=kron(A,ones(k,1));
%B=kron(ones(m,1), B);
%B=repmat(B,m,1); % no big speed improvement over 'kron', but hey...
%H=A.*B;

H = kron(A,ones(k,1)).*repmat(B,m,1);