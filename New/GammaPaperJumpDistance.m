% function GammaPaperJumpDistance
hw = 31;
deltaTime = [1 2 4 8 16 32]; % 30 ;
Color = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
figure
turnAngle = 0;
for delta = 1 % :3 % length(deltaTime)
    %     subplot(1,3,delta)
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    for i = 1 % 8:13:130 % :10 % num_files % 1:num_files
        fprintf('Loading 3DBurst.mat file %s...\n', files{i});
        R = load(files{i});
        R = get_grid_firing_centre2(R);
%         center = R.WCentroids ;
        for iBurst = 1 % :size(center,2)
            %             posCenter = center{iBurst}(:,2:3);
            
%             posCenter = cell2mat(center');
%             t = posCenter(:,1);
%             posCenter = posCenter(:,2:3);
%             
%             posCenter = posCenter/40*600; %% manually modify here according to electrdoes %%
%             posCenter = posCenter/600*2*pi;
%             posCenter = exp(1i*posCenter);
            
            % spike pattern
            posCenter = R.grid.quick.centre';
            t = find(~isnan(posCenter(:,1)));
            posCenter = posCenter(t,:);
            posCenterO = posCenter;
            posCenter = posCenter/31.5*pi; %% manually modify here according to electrdoes %%
            posCenter = exp(1i*posCenter);
            Displace1Dx = wrapToPi(posCenterO(2:end,1)/31.5*pi-posCenterO(1:end-1,1)/31.5*pi);
            Displace1Dx = Displace1Dx/(2*pi)*600;
%             Displace1Dx = posCenterO(2:end,1)-posCenterO(1:end-1,1);
%             Displace1Dx = wrapToPi(Displace1Dx/31*pi)*63/(2*pi);
            
            % 1D increment
%             R = get_grid_firing_centreYL2(R); % ,'mode','quick','win_len',bin,'jump_win',jump_win/10);
%             Displace1Dx = R.grid.quick.jump_size(:,1)/31.5*300;
%             Displace2D = R.grid.quick.jump_dist/63*600;
%             Displace = [Displace;R.grid.quick.jump_dist] ;
            
%             for stepT = 1:deltaT
%                 diffT = t(stepT+1:end)-t(1:end-stepT);
%                 ind = find(diffT==deltaT);
%                 for l = 1:length(ind)
%                     DisplaceTemp = sqrt(sum((angle(posCenter(stepT+ind(l),:)./posCenter(ind(l),:))/(2*pi)*600).^2,2));
%                     Displace = [Displace;DisplaceTemp] ;
%                 end
%             end
            
            %             DisplaceTemp = sqrt(sum((angle(posCenter(deltaT+1:end,:)./posCenter(1:end-deltaT,:))/(2*pi)*600).^2,2));
            %             Displace = [Displace;DisplaceTemp] ;
            
            
            
            % 2D step length analysis
            turnRadian = turnAngle/360*2*pi ;
%             posCenterO = center{iBurst}(:,2:3);
            pathAngle = atan2((posCenterO(deltaT+1 :end,2)-posCenterO(1: end-deltaT,2))...
                ,(posCenterO(deltaT+1 :end,1)-posCenterO(1: end-deltaT,1))) ;
            tunningIdx = abs(angdiff(pathAngle))>turnRadian ;
            tunningPts1 = find(tunningIdx==1)+1 ;
            DisplaceTemp = sqrt(sum((angle(posCenter(tunningPts1(2:end),:)./posCenter(tunningPts1(1:end-1),:))/(2*pi)*600).^2,2)); % 600
            Displace = [Displace;DisplaceTemp] ; 
            
            %             plot(center{iBurst}(:,2),center{iBurst}(:,3),'r-o')
            %             hold on
            %             plot(center{iBurst}(tunningPts1,2),center{iBurst}(tunningPts1,3),'k->')
            %             xlim([0 63])
            %             ylim([0 63])
            %             next = input('\t Next figure?');
            %             close all
        end
    end
    
    %     mean(Displace)
    Displace = Displace(Displace>0);
    numPts = 80 ;
    nEdge = logspace(log10(min(Displace)),log10(max(Displace)),numPts) ; 
    %     nEdge = logspace(log10(0.01),log10(max(Displace)),numPts) ;
    [n,nEdge] = histcounts(Displace,nEdge,'normalization','cdf');
    x = (nEdge(1:end-1)+nEdge(2:end))/2;
    h(delta) = loglog(x,1-n,'.','MarkerSize',6,'color',Color(delta,:));
    hold on;
    sigDist = Displace;
    xmin = 10 % min(sigDist)
    xmax = max(sigDist)
    pf_TP = @(x,alpha) (alpha+1)*(xmax.^(alpha+1)-xmin.^(alpha+1)).^(-1).*x.^(alpha) ;
    cdf_TP = @(x,alpha) (x.^(alpha+1)-xmin.^(alpha+1))./(xmax.^(alpha+1)-xmin.^(alpha+1)) ;
    [lambdaHat1,lambdaCI] = mle(sigDist, 'pdf',pf_TP, 'start',[-4], 'lowerbound',[-5], 'upperbound', [0])
    loglog(x,1-cdf_TP(x,lambdaHat1(1)),'r','lineWidth',2) 
    %
    hold on ;
    pd2 = fitdist(sigDist,'normal')
    y2 = cdf(pd2,x);% 2*
    loglog(x,1-y2,'k--','LineWidth',2);
    hold on
    xlim([0 5e2]) % Duration4e2
    ylim([4e-4 1])
    xlabel('Step length d(\tau) (\mum)','fontsize',8)
%     xlabel('Displacement d(\tau) (\mum)','fontsize',8)
    ylabel('{\itP} (d(\tau))','fontsize',8)
    legend({'Data','Truncated Pareto','Gaussian'},'fontsize',8)
end
% ('displacement/\eta^{0.5}') ['displacement/\eta^{', num2str(eta),'}']
% legend([h(1) h(2) h(3)],{'\tau = 2 ms','\tau = 4 ms','\tau = 8 ms',},'fontsize',8,'edgecolor','none','color','none')
% end