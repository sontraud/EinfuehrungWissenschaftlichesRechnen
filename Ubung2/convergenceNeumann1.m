tic
n=[4 9 19 39 79 99 119];
h=ones(1,length(n));
err_iter = zeros(2,length(n));
h=h./(n+1);
w=3/2;
w_array = [3/2 3];
for j = 1:length(w_array)
    w = w_array(j);
    for i=1:length(n)
        [u_h, err] =solveNeumann(n(i),w,'false');
        err_iter(j,i) = err;
        fprintf('j= %i  i = %i',j,i);
    end
end
%Überprüfe numerisch Konvergenz für gesamt Fehler von O(h^2)
loglog(h,err_iter(1,:))
hold on
loglog(h, err_iter(2,:))
loglog(h,h.^2)
legend('w=3/2','w=3','Referenz');
hold off
toc