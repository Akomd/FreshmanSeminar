% ���q�t�B���^. 
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

% �^�l������. 
x = zeros(2,sampleNum);
% �ϑ��l������. 
y = zeros(1,sampleNum);
% ����l������. 
x_y = zeros(2,sampleNum);
% �ޓx������.
lambda = zeros(1,N);

% �K�E�X���������̂��߂̐ݒ�. 
MU = ens0';
SIGMA = 0.2^2 * [1 0; 0 1];
p = 1;
obj = gmdistribution(MU, SIGMA, p);

% �����A���T���u������. 
ensemble = random(obj,N)';

% �^�l�̈�Ԗڂ����v�Z. 
x(:,1) = Pendulum_Model(x0);

% �������o�̂��߂̃C���f�b�N�X. 
ind = zeros(1,N);
for i = 1:1:N
    ind(1,i) = i;
end

k = 0;
while k < sampleNum;
    x_ensemble = Pendulum_Model(ensemble);
    y(k + 1) = Observer_Model(x(:,k + 1));
    % �ޓx�v�Z�Ɩޓx�Ɋ�Â��������o.
    inexp = ( -1 / 2 ) * (y(k + 1) - H * x_ensemble) .* ( R \ (y(k + 1) - H * x_ensemble) );
    lambda = 1 / sqrt((2*pi)^2 * det(R)) * exp( inexp );
    beta = lambda / sum(lambda);
    ensemble_ind = datasample(ind,N,'Weights',beta);
    % ���A���T���u��.
    ensemble = x_ensemble(:,ensemble_ind);
    % ����l.
    x_y(:,k + 1) = sum(ensemble,2) / N;
    % �^�l.
    if k < sampleNum - 1
        x(:,k + 2) = Pendulum_Model(x(:,k + 1));
    end
    k = k + 1;
end

index = 'test';

mkdir(index);

old = cd([pwd,'\',index]);

% �`��.
xmin = -2;
xmax = 2;
ymin = -2;
ymax = -1;
% �^�l�̕`��.
fig = figure(1);
subplot(2,1,1)
plot(l * sin(x(1,:)),-l * cos(x(1,:)),'.','Color','green');
axis([xmin xmax ymin ymax]);
axis equal;
grid on;
title('�U��q�ʒu:�^�l')
% ����l�̕`��.
subplot(2,1,2)
plot(l * sin(x_y(1,:)),-l * cos(x_y(1,:)),'.','Color','red');
axis([xmin xmax ymin ymax]);
axis equal;
grid on;
title('�U��q�ʒu:����l')

fig.PaperUnits = 'points';
fig.PaperSize = [fig.Position(3) fig.Position(4)];
fig.PaperPositionMode = 'auto';
print(fig,'pendulum','-dpdf','-r0');

% �^�l�̕`��.
fig = figure(3);
subplot(2,1,1)
plot(x(1,:),x(2,:),'.','Color','green');
axis equal;
ax = axis;
grid on;
title('(x1,x2):�^�l')
% ����l�̕`��.
subplot(2,1,2)
plot(x_y(1,:),x_y(2,:),'.','Color','red');
axis(ax);
grid on;
title('(x1,x2):����l')
fig.PaperUnits = 'points';
fig.PaperSize = [fig.Position(3) fig.Position(4)];
fig.PaperPositionMode = 'auto';
print('(x1,x2)','-dpdf','-r0');



xmax = sampleTime;
xmin = 0;
ymax = 0.35;
ymin = -0.55;
% �^�l�̕`��.
fig = figure(5);
plot((0:T:sampleTime-T),x(1,:),'.','Color','green');
axis([xmin xmax ymin ymax]);
grid on;
hold on;
% ����l�̕`��.
plot((0:T:sampleTime-T),x_y(1,:),'.','Color','red');
title('�p�x�i�^�l�F��, ����l�F�ԁj')
fig.PaperUnits = 'points';
fig.PaperSize = [fig.Position(3) fig.Position(4)];
fig.PaperPositionMode = 'auto';
print('angle','-dpdf','-r0');


xmax = sampleTime;
xmin = 0;
ymax = 0.8;
ymin = -0.8;
% �^�l�̕`��.
fig = figure(6);
plot((0:T:sampleTime-T),x(2,:),'.','Color','green');
axis([xmin xmax ymin ymax]);
grid on;
hold on;
% ����l�̕`��.
plot((0:T:sampleTime-T),x_y(2,:),'.','Color','red');
title('�p���x�i�^�l�F��, ����l�F�ԁj')
fig.PaperUnits = 'points';
fig.PaperSize = [fig.Position(3) fig.Position(4)];
fig.PaperPositionMode = 'auto';
print('vangle','-dpdf','-r0');

cd(old);

