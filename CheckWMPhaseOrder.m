% function CheckWMPhaseOrder
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
NumP = 3;
bin = 40;
Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
LoalNeu = cell(1,NumP);
p = cell(10,NumP); % num_files
R = load(files{1});
for i = 1:NumP
    dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
    LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
end
% Butterworth filter
order = 4; % 4th order
lowFreq_br = 4; % theta band
hiFreq_br = 7;
fs = 1e4;
Wn = [lowFreq_br hiFreq_br]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
figure
for no = 1:3
    for id_out = 11:20 % num_files
        fprintf('\t File name: %s\n', files{id_out});
        R = load(files{id_out});
        r = sum(movsum(full(R.spike_hist{1}(LoalNeu{no},:)),bin,2));
        LFP_theta = angle(hilbert(filter(b,a,R.LFP{1}(no,:))));
        for i = 2.26e4:length(r)/5-bin/2 % ind % 6.46e4 % length(r) %
            if r(i) >= 30
                p{id_out-10,no} = [p{id_out-10,no} LFP_theta(i)];
            end
        end
    end
    subplot(2,2,no)
    polarhistogram([p{:,no}])
end
subplot(2,2,4)
errorbar(1:3,mean(180-rad2deg(cellfun(@mean,p))),std(180-rad2deg(cellfun(@mean,p))),'o-')
hold on
errorbar(1:3,mean(180-rad2deg(cellfun(@median,p))),std(180-rad2deg(cellfun(@median,p))),'>-')
legend('mean','median')
xlabel('Loading order')
ylabel('Locked phase')
% end