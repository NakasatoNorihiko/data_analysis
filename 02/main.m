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
# pause();

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
# pause();

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
# pause();

# 交差確認法（コピペ）
clear all; rand('state', 0); randn('state', 0);
n = 50; N = 1000;
x=linspace(-3, 3, n)'; X=linspace(-3, 3, N)';
pix=pi*x; y = sin(pix)./(pix) + 0.1*x + 0.2*randn(n,1);

x2=x.^2; xx=repmat(x2,1,n)+repmat(x2',n,1)-2*x*x';
hhs=2*[0.03 0.3 3].^2; ls=[0.001 0.1 100];
# データをランダムに5分割する
m=5; u=floor(m*[0:n-1]/n)+1; u=u(randperm(n));
u
for hk=1:length(hhs)
    hh=hhs(hk); k=exp(-xx/hh);
    for i=1:m
	# 順番にデータを除外する
        ki=k(u~=i,:); kc=k(u==i,:); yi=y(u~=i); yc=y(u==i);
        for lk=1:length(ls)
            l=ls(lk); t=(ki'*ki+l*eye(n))\(ki'*yi); fc=kc*t;
            g(hk,lk,i)=mean((fc-yc.^2));
        end
    end
end
abs(mean(g,3))
[gl,ggl]=min(abs(mean(g,3)),[],2); [ghl, gghl]=min(gl);
L=ls(ggl(gghl)); HH=hhs(gghl);

K=exp(-(repmat(X.^2,1,n)+repmat(x2',N,1)-2*X*x')/HH);
k=exp(-xx/HH); t=(k^2+L*eye(n))\(k*y); F=K*t;

figure(1); clf; hold on; axis([-2.8 2.8 -1 1.5]);
plot(X, F, 'g-');
plot(x, y, 'bo');
# pause();

# 交差確認法
clear all; rand('state', 0); randn('state', 0);
n = 50; N = 1000;
x=linspace(-3, 3, n)'; X=linspace(-3, 3, N)';
pix=pi*x; y = sin(pix)./(pix) + 0.1*x + 0.2*randn(n,1); truey = sin(pix)./(pix) + 0.1*x

x2=x.^2; xx=repmat(x2,1,n)+repmat(x2',n,1)-2*x*x';
hhlist=2*[0.001 0.01 0.1 1 10].^2; llist=[0.001 0.01 0.1 1 10];
# データをランダムに5分割する
m=5; u=floor(m*[0:n-1]/n)+1; u=u(randperm(n));
for hk=1:length(hhlist)
    hh=hhlist(hk); k=exp(-xx/hh);
    for i=1:m
	# 順番にデータを除外する
        ki=k(u~=i,:); kc=k(u==i,:); yi=y(u~=i); yc=y(u==i);
        for lk=1:length(llist)
            l=llist(lk); t=(ki'*ki+l*eye(n))\(ki'*yi); fc=kc*t;
            g(hk,lk,i)=mean((fc-yc).^2);
        end
    end
end
mean(g,3)
[gl,ggl]=min(mean(g,3),[],2); [ghl, gghl]=min(gl);
L=llist(ggl(gghl)); HH=hhlist(gghl);

K=exp(-(repmat(X.^2,1,n)+repmat(x2',N,1)-2*X*x')/HH);
k=exp(-xx/HH); t=(k^2+L*eye(n))\(k*y); F=K*t;

figure(1); clf; hold on; axis([-2.8 2.8 -1 1.5]);
plot(X, F, 'g-');
plot(x, y, 'bo');
plot(x, truey, 'r-');
print("plot.png");
