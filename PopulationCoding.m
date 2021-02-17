% function PopulationCoding
figure_width = 11.4; % cm
figure_hight = 11.4; % cm
figure('NumberTitle','off','name', 'WMPopulationCoding', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
bin = 500; % 4ms
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
xr = [];
yr = [];
no = 1;
Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
LoalNeu = cell(1,3);
R = load(files{1});
for i = 1:3
    dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
    LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
end
for id_out = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    StiNeu = LoalNeu;        
    LoadC = [];
    RecallC = [];
    r = sum(movsum(full(R.spike_hist{1}(StiNeu{no},:)),bin,2));
    candx = Lattice(StiNeu{no},1);
    candy = Lattice(StiNeu{no},2);
    for i = 2.26e4:length(r)-bin/2 % 4.46e4 % length(r) % 
        if r(i) >= 120 % >25
            % neurons within WM area
            spike_x_pos_o = repmat(candx,1,bin).*R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2);
            spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2));
            spike_y_pos_o = repmat(candy,1,bin).*R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2);
            spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2));
            [x,y,width,~,height,~] = fit_bayesian_bump_2_spikes_circular(spike_x_pos_o,spike_y_pos_o,2*hw+1,'quick');
            subplot(2,2,1)
            plot(spike_x_pos_o,spike_y_pos_o,'b.')
            xlabel('x','fontSize',10)
            ylabel('y','fontSize',10)
            text(-0.3,1,'A','Units', 'Normalized','FontSize',12)
            % xlim([-0.6,0.6]) % ([-1,1]) % ([-17 -15])
            % ylim([-0.6,0.6]) % ([-1,1]) % ([-17 -15])
            subplot(2,2,2)
            edges = floor(min(spike_y_pos_o)):ceil(max(spike_y_pos_o)); % y-3*width:y+3*width;
            [N,edges] = histcounts(spike_y_pos_o,edges,'Normalization','probability'); % ,edges
            edges = (edges(1:end-1)+edges(2:end))/2;
            plot(edges,N,'o')
            hold on
            ycos = -0.01*height*cos(2*pi/(6*width)*(edges+y))+0.005*height;
            plot(edges,ycos,'-','LineWidth',1.5);
            xlabel('y Coordinate','fontSize',10)
            ylabel('Cosine fit','fontSize',10)
            text(-0.3,1,'B','Units', 'Normalized','FontSize',12)
            subplot(2,2,3)
            edges = floor(min(spike_x_pos_o)):ceil(max(spike_x_pos_o)); % x-3*width:x+3*width;
            [N,edges] = histcounts(spike_x_pos_o,edges,'Normalization','probability'); % ,edges
            edges = (edges(1:end-1)+edges(2:end))/2;
            plot(edges,N,'o')
            hold on
            xcos = -0.01*height*cos(2*pi/(6*width)*(edges+x))+0.005*height;
            plot(edges,xcos,'-','LineWidth',1.5);
            xlabel('x Coordinate','fontSize',10)
            ylabel('Cosine fit','fontSize',10)
            text(-0.3,1,'C','Units', 'Normalized','FontSize',12)
            next = input('\t Next figure?');
            close all
        end
    end
end

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc WMPopulationCoding
% end