tic
n_array=[4 9 19 39 79 99];
err_array = zeros(1,length(n_array));
w=1;
h_array=ones(1,length(n_array));
h_array=h_array./(n_array+1);
%wir brauchen nur A_h+ und R
for i = 1:length(n_array)
    h=h_array(i);
    n=n_array(i);
    T=sparse(1:n+2,1:n+2,4*ones(1,n+2))+sparse(2:n+2,1:n+1,-ones(1,n+1),n+2,n+2)+sparse(1:n+1,2:n+2,-ones(1,n+1),n+2,n+2);
    T(1,2)=-2;
    T(n+2,n+1) = -2;
    T_schlange = sparse(2:n+2,1:n+1,-ones(1,n+1),n+2,n+2)+sparse(1:n+1,2:n+2,-ones(1,n+1),n+2,n+2);
    T_schlange(1,2) = -2;
    T_schlange(n+2,n+1) = -2;
    %kron_1 = kron(speye(n+2),T);
    %kron_2 = kron(T_schlange,speye(n+2));
    %Erstelle A_h
    A = 1/h^2*((kron(speye(n+2),T))+kron(T_schlange,speye(n+2)));
    clear T,T_schlange;
    %A = 1/h^2*(kron_1+kron_2);
    r=zeros(n+2,n+2);
    r(1,:) = w;
    r(:,n+2) = w;
    r(n+2,:) = w;
    r(:,1) = w;
    %Erstelle R
    R=sparse(1:(n+2)^2,1:(n+2)^2,r,(n+2)^2,(n+2)^2);
    %A_fehler = pinv(full(A))*full(R);
    err_array(i)=norm(pinv(full(A))*full(R),'inf');
    clear A,R;
    fprintf('i=%i mit Fehler %f \n',i,err_array(i));
end
%Überprüfe numerisch Konvergenz für ||A_h+R|| von O(h)
loglog(h_array,err_array)
hold on
loglog(h_array,h_array)
loglog(h_array,h_array.*h_array);
legend('w=1','Referenz O(h)', 'Referenz O(h^2)');
toc

