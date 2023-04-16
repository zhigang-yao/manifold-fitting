function T = Res2Tab(Times1,Times2)
[m,n] = size(Times1);
T1 = zeros(m*n,4);
t = 0;
for ii = 1:m
    for jj = 1:n
        t = t + 1;
        T1(t,:) = [1,ii,jj,Times1(ii,jj)];
    end
end

[m,n] = size(Times2);
T2 = zeros(m*n,4);
t = 0;
for ii = 1:m
    for jj = 1:n
        t = t + 1;
        T2(t,:) = [2,ii,jj,Times2(ii,jj)];
    end
end

T = [T1;T2];

T = array2table(T, "VariableNames",["Method","NumSample","Rep",'Value']);
