% 粒子フィルタ. 
global g;
global l;
global m;
global T;
global Q;
global R;
global H;
g = 9.8;
l = 2.0;
m = 1.0;
T = 0.001;
Q = 0.2^2;
R = 0.1^2;
H = [0 1];

sampleTime = 10;
sampleNum = int32( sampleTime / T );

x0 = [pi/24;0];
ens0 = [pi/24; 0];
N = 100;

% 真値初期化. 
x = zeros(2,sampleNum);
% 観測値初期化. 
y = zeros(1,sampleNum);
% 推定値初期化. 
x_y = zeros(2,sampleNum);
% 尤度初期化.
lambda = zeros(1,N);

% ガウス乱数生成のための設定. 
MU = ens0';
SIGMA = 0.2^2 * [1 0; 0 1];
p = 1;
obj = gmdistribution(MU, SIGMA, p);

% 初期アンサンブル生成. 
ensemble = random(obj,N)';

% 真値の一番目だけ計算. 
x(:,1) = Pendulum_Model(x0);

% 復元抽出のためのインデックス. 
ind = zeros(1,N);
for i = 1:1:N
    ind(1,i) = i;
end

k = 0;
while k < sampleNum;
    x_ensemble = Pendulum_Model(ensemble);
    y(k + 1) = Observer_Model(x(:,k + 1));
    % 尤度計算と尤度に基づく復元抽出.
    inexp = ( -1 / 2 ) * (y(k + 1) - H * x_ensemble) .* ( R \ (y(k + 1) - H * x_ensemble) );
    lambda = 1 / sqrt((2*pi)^2 * det(R)) * exp( inexp );
    beta = lambda / sum(lambda);
    ensemble_ind = datasample(ind,N,'Weights',beta);
    % 次アンサンブル.
    ensemble = x_ensemble(:,ensemble_ind);
    % 推定値.
    x_y(:,k + 1) = sum(ensemble,2) / N;
    % 真値.
    if k < sampleNum - 1
        x(:,k + 2) = Pendulum_Model(x(:,k + 1));
    end
    k = k + 1;
end

index = 'test';

mkdir(index);

old = cd([pwd,'\',index]);

% 描画.
xmin = -2;
xmax = 2;
ymin = -2;
ymax = -1;
% 真値の描画.
fig = figure(1);
subplot(2,1,1)
plot(l * sin(x(1,:)),-l * cos(x(1,:)),'.','Color','green');
axis([xmin xmax ymin ymax]);
axis equal;
grid on;
title('振り子位置:真値')
% 推定値の描画.
subplot(2,1,2)
plot(l * sin(x_y(1,:)),-l * cos(x_y(1,:)),'.','Color','red');
axis([xmin xmax ymin ymax]);
axis equal;
grid on;
title('振り子位置:推定値')

fig.PaperUnits = 'points';
fig.PaperSize = [fig.Position(3) fig.Position(4)];
fig.PaperPositionMode = 'auto';
print(fig,'pendulum','-dpdf','-r0');

% 真値の描画.
fig = figure(3);
subplot(2,1,1)
plot(x(1,:),x(2,:),'.','Color','green');
axis equal;
ax = axis;
grid on;
title('(x1,x2):真値')
% 推定値の描画.
subplot(2,1,2)
plot(x_y(1,:),x_y(2,:),'.','Color','red');
axis(ax);
grid on;
title('(x1,x2):推定値')
fig.PaperUnits = 'points';
fig.PaperSize = [fig.Position(3) fig.Position(4)];
fig.PaperPositionMode = 'auto';
print('(x1,x2)','-dpdf','-r0');



xmax = sampleTime;
xmin = 0;
ymax = 0.35;
ymin = -0.55;
% 真値の描画.
fig = figure(5);
plot((0:T:sampleTime-T),x(1,:),'.','Color','green');
axis([xmin xmax ymin ymax]);
grid on;
hold on;
% 推定値の描画.
plot((0:T:sampleTime-T),x_y(1,:),'.','Color','red');
title('角度（真値：緑, 推定値：赤）')
fig.PaperUnits = 'points';
fig.PaperSize = [fig.Position(3) fig.Position(4)];
fig.PaperPositionMode = 'auto';
print('angle','-dpdf','-r0');


xmax = sampleTime;
xmin = 0;
ymax = 0.8;
ymin = -0.8;
% 真値の描画.
fig = figure(6);
plot((0:T:sampleTime-T),x(2,:),'.','Color','green');
axis([xmin xmax ymin ymax]);
grid on;
hold on;
% 推定値の描画.
plot((0:T:sampleTime-T),x_y(2,:),'.','Color','red');
title('角速度（真値：緑, 推定値：赤）')
fig.PaperUnits = 'points';
fig.PaperSize = [fig.Position(3) fig.Position(4)];
fig.PaperPositionMode = 'auto';
print('vangle','-dpdf','-r0');

cd(old);

