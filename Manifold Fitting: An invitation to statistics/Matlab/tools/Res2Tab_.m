function T = Res2Tab_(Times1)
[m,n] = size(Times1);
T1 = zeros(m*n,3);
t = 0;
for ii = 1:m
    for jj = 1:n
        t = t + 1;
        T1(t,:) = [ii,jj,Times1(ii,jj)];
    end
end

T = T1;

T = array2table(T, "VariableNames",["NumSample","Rep",'Value']);
