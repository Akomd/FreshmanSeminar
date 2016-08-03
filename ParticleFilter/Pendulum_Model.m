% U‚èqƒ‚ƒfƒ‹. 
function y = Pendulum_Model(x)
global g;
global l;
global m;
global T;
global Q;
sigma = sqrt(Q);
v = randn(1,size(x,2)) * sigma;
y = x + ( [x(2,:); -g / l * sin(  x(1,:) )] + [0; -1 / ( m * l )] * v ) * T;
end
