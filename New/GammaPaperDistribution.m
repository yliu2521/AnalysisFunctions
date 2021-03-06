% function GammaPaperDistribution
sigDist = Duration(Duration>59); % Interval
xmin = min(sigDist) % 95* 18*
xmax = max(sigDist) 
pf_TP = @(x,alpha) (alpha+1)*(xmax.^(alpha+1)-xmin.^(alpha+1)).^(-1).*x.^(alpha) ;
cdf_TP = @(x,alpha) (x.^(alpha+1)-xmin.^(alpha+1))./(xmax.^(alpha+1)-xmin.^(alpha+1)) ;
[lambdaHat1,lambdaCI] = mle(sigDist, 'pdf',pf_TP, 'start',[-4], 'lowerbound',[-5], 'upperbound', [0])
figure;
numPts = 200 ;
nEdge = logspace(log10(min(sigDist)),log10(max(sigDist)),numPts) ; 
[x,n] = histcounts(sigDist,nEdge,'normalization','cdf') ;
n0 = 0.5*((n(2:end))+(n(1:end-1))) ;
loglog(n0,1-x,'b.','MarkerSize',12)
hold on
n1 = n0(n0>50);
loglog(n0,1-cdf_TP(n0,lambdaHat1(1)),'r','lineWidth',2) 
pd2 = fitdist(sigDist','normal')
y2 = cdf(pd2,n0);% 2*
hold on
loglog(n0,1-y2,'k--','LineWidth',2);
xlim([0 4e2]) % Duration4e2 Interval7e2
ylim([6e-5 1])
% set(gca,'XScale', 'log'); 
legend({'Data','Truncated Pareto','Gaussian'},'fontsize',8)
xlabel('Duration','fontsize',8)
ylabel('Probability {\itP} (s{\itT})','fontsize',8)
% end