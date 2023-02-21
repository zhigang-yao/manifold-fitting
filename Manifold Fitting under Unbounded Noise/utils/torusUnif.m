function X = torusUnif(n, R, r)
% generate n samples uniformly from a torus
% R : the radius from the center of the hole to the center of the torus tube
% r : the radius of the torus tube. 
%
% There is an R package with the same name
% Algorithm 1 of Diaconis P, Holmes S, and Shahshahani M (2013). "Sampling from a manifold." 
% Advances in Modern Statistical Theory and Applications: A Festschrift in honor of Morris L. Eaton. 
% Institute of Mathematical Statistics, 102-125. 

count = 0;
theta = -ones(1,n);

while count < n
    xvec = rand(1)*2*pi;
    yvec = rand(1)/pi;
    fx = (1+r/R*cos(xvec)) / (2*pi);
    if yvec < fx
        count = count + 1;
        theta(count) = xvec;
    end
end

phi = rand(1,n)*2*pi;
x = (R + r*cos(theta)) .* cos(phi);
y = (R + r*cos(theta)) .* sin(phi);
z = r*sin(theta);

X = [x;y;z];
