% ŠÏ‘ªƒ‚ƒfƒ‹. 
function y = Observer_Model(x)
global R;
global H;
sigma = sqrt(R);
w = randn(1,size(x,2)) * sigma;
y = H * x + w;
end