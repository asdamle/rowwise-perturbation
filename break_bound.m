clear 
close all
clc

% try some stuff for matrix perturbations to break our bound

setN = 2*round(logspace(4,6,20)/2);
numN = length(setN);


for k = 1:numN
    n = setN(k)
    v = ones(n,1);
    v = v/norm(v);
    c = 1/21;
    Et = sparse(n,n);
    w = ones(n,1);
    w = w/norm(w);
    Et(1,:) = w';
    Et = c*Et;
    yB = ones(n,1);
    yB(2:end) = ((sqrt(n)-2*c)/(c*(n-1)))*yB(2:end);
    idx = round(1*n/2);
    c*w'*yB
    yB(idx:end) = -1*yB(idx:end);
    Et(1,idx:end) = -1*Et(1,idx:end);
    Et = Et+Et';    
    z  = yB - Et*yB;
    
    V2x = @(x) x - v*(v'*x);
    E22 = @(x) V2x(Et*V2x(x));
    gamma = (n^(1/3));
    Ex = @(x) V2x(z)*(v'*x)/gamma + v*(V2x(z)'*x)/gamma + (Et*x + v*(v'*Et*v)*(v'*x) - v*((v'*Et)*x) - (Et*v)*(v'*x));
    
    Ex_svd = @(x) [Ex(x((n+1):2*n));Ex(x(1:n))];
    tt(k) = max(abs(eigs(@(x) Ex_svd(x),2*n,4)));
    
    Exn = @(x) (Ex(x)/(tt(k))/(n^(1/3)));
    
    Exn_svd = @(x) [Exn(x((n+1):2*n));Exn(x(1:n))];
    eN(k) = max(abs(eigs(@(x) Exn_svd(x),2*n,4)));
    
    [VV, DD] = eigs(@(x) v*(v'*x) + Exn(x),n,1);
    
    vt = VV(:,1);
    yt = vt - v*(v'*vt);
    eS(k) = sqrt(1-(v'*vt)^2);
    eI(k) = min(norm(v-vt,'inf'),norm(v+vt,'inf'));
    e21(k) = norm(V2x(Exn(v)),'inf');
    b1 = zeros(n,1); b1(1) = 1;
    e2Inf(k) = norm(V2x(Exn(b1)));
    eVEY(k) = norm(V2x(Exn(yt)),'inf');
    eEY_split(k) = e2Inf(k)*norm(yt);
    eEE(k) = e2Inf(k)*eN(k);
    e1(k) = norm(v,'inf')*eN(k)^2;
    
    usign = sign(v'*vt);
    eV1(k) = norm(v - usign*v*(vt'*v),'inf');
    eV2(k) = norm((v-usign*vt)-v*(v'*(v-usign*vt)),'inf');
    
    
end

figure;
a = loglog(setN,eV1,'-','LineWidth',5);
set(a,'Marker','.','MarkerSize',40)
hold on
loglog(setN,e1,'--','LineWidth',5,'MarkerSize',10);
legend({'$\|V_1V_1^T(V_1-\hat{V}_1\tilde{U})\|_{2,\infty}$','$\|V_1\|_{2,\infty}\|E\|_2^2$'},'Interpreter', 'Latex');
xlim([setN(1) setN(end)])
xlabel('n','Interpreter','Latex')
set(gca,'FontSize',58)
set(gcf,'position',[0 0 1200 1200])

figure;
a = loglog(setN,eV2,'-','LineWidth',5);
set(a,'Marker','.','MarkerSize',40)
hold on
loglog(setN,e21,'-.','LineWidth',5,'MarkerSize',10);
loglog(setN,eEE,'-.','LineWidth',5,'MarkerSize',10);
legend({'$\|V_2V_2^T(V_1-\hat{V}_1\tilde{U})\|_{2,\infty}$','$\|V_{2}E_{21}\|_{2,\infty}$','$\|V_2V_2^TE\|_{2,\infty}\|E\|_2$'},'Interpreter', 'Latex');
xlim([setN(1) setN(end)])
xlabel('n','Interpreter','Latex')
ylim([10^-5/2,10^-3])
set(gca,'FontSize',58)
set(gcf,'position',[0 0 1200 1200])