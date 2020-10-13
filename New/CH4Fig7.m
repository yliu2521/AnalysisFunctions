% function CH4Fig7
figure_width = 11.4; % cm
figure_hight = 11.4; % cm
figure('NumberTitle','off','name', 'CH4Fig7', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
dir_strut2 = dir('*_config_data.mat');
num_files2 = length(dir_strut2);
files2 = cell(1,num_files2);
for id_out = 1:num_files2
    files2{id_out} = dir_strut2(id_out).name;
end
bin = 40; % 4ms
FR = zeros(1,num_files);
Duration = zeros(1,num_files);
for i = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', i, num_files);
    fprintf('\t File name: %s\n', files{i});
    R = load(files{i});
    load(files2{i},'StiNeu')
%     FR(i) = mean(R.Analysis.rate{1});
%     FR(i) = sum(sum(R.spike_hist{1}(StiNeu{1},2.25e4+1:end)))/50/17.75; % Hz
    FR(i) = sum(R.num_spikes{1})/3969/(length(R.num_spikes{1}/1e4)); % 20s or 10s
    r = sum(movsum(full(R.spike_hist{1}(StiNeu{1},:)),bin,2));
    Dind = find(r > 25);
    if Dind(end) > 2.25e4
        Duration(i) = Dind(end)*0.1-2.25e3;
    end
end
FR = vec2mat(FR,9);
fr = mean(FR);
frSTD = std(FR);
Fr = 0.6:0.1:1.4;
v1 = polyfit(log10(Fr(1:5)),log10(fr(1:5)),1);
x1 = Fr(1:7);
y1 = 10^v1(2)*x1.^v1(1);
v2 = polyfit(log10(Fr(5:end)),log10(fr(5:end)),1);
x2 = Fr(3:end);
y2 = 10^v2(2)*x2.^v2(1);
subplot(1,2,1)
errorbar(Fr,fr,frSTD,'o') % ,'MarkerSize',6,'CapSize',6,'LineWidth',1.5)
set(gca,'YScale','log')
hold on
semilogy(x1,y1,'LineWidth',1.5)
hold on
semilogy(x2,y2,'LineWidth',1.5)
xlabel('Scaled IE ratio','fontsize',10)
ylabel('Firing Rate(Hz)','fontsize',10)
text(-0.1,1,'A','Units', 'Normalized','FontSize',12)
axes('Position',[.3 .6 .15 .15])
box on
errorbar(Fr,fr,frSTD,'o')
hold on
plot(x1,y1,'LineWidth',1.5)
hold on
plot(x2,y2,'LineWidth',1.5)
% -v1(1)
% -v2(1)
% xlabel('Scaled IE ratio','fontsize',10)
% ylabel('Firing Rate(Hz)','fontsize',10)

Duration = Duration*1e-3; % s
Duration = vec2mat(Duration,9);
dur = mean(Duration);
durSTD = std(Duration);
v2 = polyfit(log10(Fr(5:end)),log10(dur(5:end)),1);
x2 = Fr(4:end);
y2 = 10^v2(2)*x2.^v2(1);
subplot(1,2,2)
errorbar(Fr,dur,durSTD,'o-.')
xlabel('Scaled IE ratio','fontsize',10)
ylabel('Working Memory Duration(s)','fontsize',10)
text(-0.1,1,'B','Units', 'Normalized','FontSize',12)

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH4Fig7 % this is the trick!!
% end