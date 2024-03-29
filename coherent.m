clear 
close all
clc
rng(1)
% try some stuff for matrix perturbations
disp('coherent test')


setN = 6*round(logspace(3,4,10)/4);
numN = length(setN);
trials = 30;


eS = zeros(trials,numN);
eP = zeros(trials,numN);
e1 = zeros(trials,numN);
e2 = zeros(trials,numN);

s = 1

f = @(x) 1./(x.^(s));

for l = 1:numN
    n = setN(l)
    k = 2;
    U = ones(n,2);
    U(1:n/2,2) = -1;
    U = U * diag(1./sqrt(sum(U.^2)));
    z = zeros(n,1);
    z(1) = 1;
    z(2) = -1;
    z = z/norm(z);
    for p = 1:trials
        E = triu(f(n)*randn(n,n));
        E = E + triu(E,1)';
        A = (4*U)*U'+z*z'+E;
        [V,D] = eigs(A,2);

        [Ut, St, Vt] = svd(V'*U);
        W = Ut*Vt';
        V2x = @(x) x - U*(U'*x);
        eS(p,l) = norm(U-V*W,'fro');
        eP(p,l) = norm2inf(U - V*W);

    end
end



%plotting
M_eS = mean(eS);
M_eP = mean(eP);
M_e1 = mean(e1);
M_e2 = mean(e2);

qlow = 0.05;
qhigh = 0.95;

MSl = quantile(eS,qlow);
MSh = quantile(eS,qhigh);

MPl = quantile(eP,qlow);
MPh = quantile(eP,qhigh);


gs = @(x) sqrt(x).*f(x);
g = @(x) sqrt(log(x)).*(f(x));

Gn = g(setN);
Gsn = gs(setN);

Gn = (1/Gn(1))*(1.2*M_eP(1))*Gn;
Gsn = (1/Gsn(1))*(1.2*M_eS(1))*Gsn;

figure;
loglog(setN,M_eS,'-*','LineWidth',5,'MarkerSize',10);
hold on
loglog(setN,M_eP,'-*','LineWidth',5,'MarkerSize',10);
loglog(setN,Gsn,'-.','LineWidth',5);
loglog(setN,Gn,'-.','LineWidth',5);


alpha = 0.15;
C = get(gca,'ColorOrder');
fill([setN fliplr(setN)],[MSl fliplr(MSh)],C(1,:),'FaceAlpha',alpha,'LineStyle','none')
fill([setN fliplr(setN)],[MPl fliplr(MPh)],C(2,:),'FaceAlpha',alpha,'LineStyle','none')

if s == 1
    legend({'$\|V-\widehat{V}\widetilde{U}\|_{F}$','$\|V-\widehat{V}\widetilde{U}\|_{2,\infty}$','$1/n^{1/2}$','$\sqrt{\log n}/n$'},'Interpreter', 'Latex')
elseif s == 3/4
    legend({'$\|V-\widehat{V}\widetilde{U}\|_{F}$','$\|V-\widehat{V}\widetilde{U}\|_{2,\infty}$','$1/n^{1/4}$','$\sqrt{\log n}/n^{3/4}$'},'Interpreter', 'Latex')
else
    disp('need legend')
end
xlim([setN(1) setN(end)])
xlabel('n','Interpreter','Latex')
set(gca,'FontSize',58)
set(gcf,'position',[0 0 1200 1200])




function x = norm2inf(A)
    x = max(sqrt(sum(A.^2,2)));
end






