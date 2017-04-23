# 三角多項式モデル
clear all; rand('state', 0); randn('state', 0);
n=50; N=1000;
x=linspace(-3, 3, n)'; X=linspace(-3, 3, N)';
pix=pi*x; y=sin(pix)./(pix)+0.1*x+0.05*randn(n,1);

p(:,1)=ones(n,1); P(:,1)=ones(N,1);
for j=1:15
    p(:, 2*j)=sin(j/2*x); p(:, 2*j+1)=cos(j/2*x);
    P(:, 2*j)=sin(j/2*X); P(:, 2*j+1)=cos(j/2*X);
end
t=(p'*p)\(p'*y); F=P*t;

figure(1); clf; hold on; axis([-2.8 2.8 -0.5 1.2]);
plot(X, F, 'g-'); plot(x, y, 'bo');
pause();

# （ガウス）カーネルモデル 通常の最小二乗回帰 
clear all; rand('state', 0); randn('state', 0);
n = 50; N = 1000;
x=linspace(-3, 3, n)'; X = linspace(-3, 3, N)';
pix = pi * x; y = sin(pix)./(pix) + 0.1*x + 0.2*randn(n,1);

x2=x.^2; X2=X.^2; hh = 2*0.3^2; l = 0.1;
k=exp(-(repmat(x2,1,n)+repmat(x2',n,1)-2*x*x')/hh);
K=exp(-(repmat(X2,1,n)+repmat(x2',N,1)-2*X*x')/hh);
t=(k)\(y);
F=K*t;

figure(1); clf; hold on; axis([-2.8 2.8 -1 1.5]);
plot(X, F, 'g-');
plot(x, y, 'bo');
pause();

# （ガウス）カーネルモデル l2-制約付き最小二乗回帰 
clear all; rand('state', 0); randn('state', 0);
n = 50; N = 1000;
x=linspace(-3, 3, n)'; X = linspace(-3, 3, N)';
pix = pi * x; y = sin(pix)./(pix) + 0.1*x + 0.2*randn(n,1);

x2=x.^2; X2=X.^2; hh = 2*0.3^2; l = 0.1;
k=exp(-(repmat(x2,1,n)+repmat(x2',n,1)-2*x*x')/hh);
K=exp(-(repmat(X2,1,n)+repmat(x2',N,1)-2*X*x')/hh);
t=(k^2+l*eye(n))\(k*y);
F=K*t;

figure(1); clf; hold on; axis([-2.8 2.8 -1 1.5]);
plot(X, F, 'g-');
plot(x, y, 'bo');
pause();

