clear 
close all
clc

% try some stuff for matrix perturbations to break our bound

setN = 25000 + 25000*(0:19);
% setN = [setN setN(end)+setN];
numN = length(setN);


for k = 1:numN
    n = setN(k)
    v = ones(n,1);
    v = v/norm(v);
%     [Q, R] = qr([v randn(n,n-1)]);
%     V2 = Q(:,2:n);
    c = 1/7;
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
%     z = yB;
    z  = yB - Et*yB;
    
    V2x = @(x) x - v*(v'*x);
    E22 = @(x) V2x(Et*V2x(x));
%     E22 = V2'*(Et)*V2;
%     E21 = V2'*z;
%     Eeig = [0 E21'; E21 E22];
%     E = [v V2]*Eeig*([v V2]');
%     E = V2x(z)*v' + (V2x(z)*v')' + (Et + v*(v'*Et*v)*v' - v*(v'*Et) - (Et*v)*v');
    Ex = @(x) V2x(z)*(v'*x)/2.5 + v*(V2x(z)'*x)/2.5 + (Et*x + v*(v'*Et*v)*(v'*x) - v*((v'*Et)*x) - (Et*v)*(v'*x));
    tt = sqrt(eigs(@(x) Ex(Ex(x)),n,1));
    Exn = @(x) (Ex(x)/(tt)/(n^(1/4)));
    eN(k) = sqrt(eigs(@(x) Exn(Exn(x)),n,1));
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
    e1(k) = norm(v,'inf')*eN(k);
    
    
end

figure;
loglog(setN,eI,'-*','LineWidth',5,'MarkerSize',10);
hold on
loglog(setN,e1,'--','LineWidth',5,'MarkerSize',10);
loglog(setN,e21,'-.','LineWidth',5,'MarkerSize',10);
loglog(setN,eEE,'-.','LineWidth',5,'MarkerSize',10);
legend({'$\|V_1-\hat{V}_1\tilde{U}\|_{2,\infty}$','$\|V_1\|_{2,\infty}\|E\|_2$','$\|V_{2}E_{21}\|_{2,\infty}$','$\|V_2V_2^TE\|_{2,\infty}\|E\|_2$'},'Interpreter', 'Latex');
xlim([setN(1) setN(end)])
xlabel('n','Interpreter','Latex')
ylim([5*10^-5 .7*10^-3])
set(gca,'FontSize',56)