function DistributedCommunication1(varargin)
% select electrode No as send and others as receive;
% work on bins which MUST contain patterns.
tic;

%%% read all RYG.mat %%%
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end

dir_strut2 = dir('UPattern-0*.mat');
% dir_strut2 = dir('AmpPatternLFPCut-0*.mat');
num_files2 = length(dir_strut2);
files2 = cell(1,num_files2);
for id_out = 1:num_files2
    files2{id_out} = dir_strut2(id_out).name;
end

%%% basic setting %%%
Opt.taux = -1; % arbitrary value>= tauy
Opt.trperm = 0; % number of trials - 1
Opt.method = 'dr';
Opt.bias = 'qe';
Opt.nt = 1; % number of calculating trails
Opt.testMode = 1;
[Lattice2, ~] = lattice_nD(2, 9.5); % creat 4*4

%%% Loop number for PBS array job %%%
loop_num = 0;

bin = 300; % 30ms
nonoverlap = bin; % unit: 0.1ms
start = 10006; % R.grid.t_mid = 26:10:99966
No = 1; % No. of Electrode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:num_files
    
    % For PBS array job
    loop_num = loop_num + 1;
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    
    fprintf('Loading RYG.mat file %s...\n', files{i});
    R = load(files{i});
    fprintf('Loading RYG.mat file %s...\n', files2{i-7});
    P = load(files2{i-7});
    [no,steps] = size(R.LFP.LFP_broad);
    LFP_certain = zeros(no,steps);
%     LFP_suffle = zeros(no,steps);
    pt = 1e4*P.ts; % P.ts: second  pt: time step
    dc = 0;
    
    % set seed
    date_now = datestr(now,'yyyymmddHHMM-SSFFF');
    scan_temp = textscan(date_now,'%s','Delimiter','-');
    rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
    alg='twister';
    rng(rand_seed,alg); %  Note that this effect is global!
    
    fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
    
    % Butterworth filter
    order = 4; % 4th order
    lowFreq_br = 5;
    hiFreq_br = 100;
    Wn = [lowFreq_br hiFreq_br]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
    for n = 1:no
        LFP_certain(n,:) = filter(b,a,R.LFP.LFP{1}(n,:)); % R.LFP.LFP{1}
%         LFP_suffle(n,:) = ifft(fft(LFP_certain(n,:)).*exp(1j*pi*rand(1,steps)),'symmetric');
    end
%     baseline = min([LFP_certain(:);LFP_suffle(:)]);
%     baseline = abs(floor(baseline));
%     LFP_certain = LFP_certain + baseline;
%     LFP_suffle = LFP_suffle + baseline;
    Hilbert = zeros(size(LFP_certain));
    for n = 1:no
        Hilbert(n,:) = hilbert(LFP_certain(n,:));
    end
%     send = No;
%     d = lattice_nD_find_dist(Lattice2,9.5,send);
%     receive = find(d > 7)'; % d > 1.5 % d > 0 
%     receive = datasample(receive,4,'Replace',false);
    for tt = start:nonoverlap:steps-2*bin
        period = tt:(tt+bin-1);
        ind = find(ismember(period,pt));
        if ~isempty(ind)
            t = period(ind(1));
            ind = find(pt == t);
            d = lattice_nD_find_dist(Lattice2,9.5,P.x(ind),P.y(ind));
            send = find(d == min(d));
            send = send(1);
            receive = find(d > P.r(ind));
            receive = datasample(receive,7,'Replace',false);
        else
            continue
        end
        for k = 1:length(receive)
            nte = zeros(1,bin);
            ntesh = zeros(1,bin);
            % LFP phase as candidate signal
            X = repmat(angle(Hilbert(send,t:(t+bin-1)))+pi,4,1);
            Y = repmat(angle(Hilbert(receive(k),t:(t+bin-1)))+pi,4,1);
            x = X(1,:);
            y = Y(1,:);
            Xsh = repmat(x(randperm(bin)),4,1);
            Ysh = repmat(y(randperm(bin)),4,1);
            % LFP as candidate signal
%             X = repmat(LFP_certain(send,t:(t+bin-1)),4,1);
%             Y = repmat(LFP_certain(receive(k),t:(t+bin-1)),4,1);
%             x = X(1,:);
%             y = Y(1,:);
%             Xsh = repmat(LFP_suffle(send,t:(t+bin-1)),4,1);
%             Ysh = repmat(LFP_suffle(receive(k),t:(t+bin-1)),4,1);
            for tau = -(1:bin)
                Opt.tauy = tau;
                try
                    [NTE] = transferentropy(X',Y',Opt,'NTE');
                    NTE = NTE(NTE ~= Inf);
                    NTE = NTE(NTE ~= -Inf);
                    [NTEsh] = transferentropy(Xsh',Ysh',Opt,'TE','NTE');
                    NTEsh = NTEsh(NTEsh ~= Inf);
                    NTEsh = NTEsh(NTEsh ~= -Inf);
                    ntesh(abs(tau)) = nanmean(NTEsh);
                    if nanmean(NTE) > 0
                        nte(abs(tau)) = nanmean(NTE);
                    end
                catch
                end
            end
            if max(nte) > max([max(ntesh),0])
                dc = dc + 1;
            end
        end
    end
    save(['NTE2Random4E-',sprintf('%04g',loop_num),'.mat'],'dc')
end
toc;
end