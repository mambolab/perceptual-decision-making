% This function calculates the time evolution of the alpha power modulation 
% in the easy and hard coherence levels locked to response onset (Figure 4B).
% 
% To make it running without any error, you need to have fieldtrip
% installed and data downloaded.
%
% Author: Alessio Basti and Roberto Guidotti
% email: roberto [dot] guidotti [at] unich [dot] it
% Version: May 2020
% Copyright (C) 2020-2022 - Methods and Models for Brain Oscillations
% (MAMBO) group, Dept. of Neuroscience, Imaging and Clinical Sciences, 
% G. d'Annunzio University Chieti-Pescara 

restoredefaultpath;

clear
clc
close all

ft_dir = '/home/mambolab/git/fieldtrip/';
addpath(ft_dir);
ft_defaults;

pn = '/home/mambolab/git/perceptual-decision-making/data/'
n_subjects = 16;
band = 'alpha';


%% Alpha power response-locked

for k = 1:n_subjects
      
    fname = ['/sub-', sprintf('%02d', k), ...
             '_task-pdm', ...
             '_band-', band, ...
             '_desc-', sprintf('%s', desc{f}), ...
             '_powspec.mat'];

    
    data = load(fullfile(pn, fname));
        
    time = data.tfa.time;
    indx_all = (time >= time(nearest(time,-1.48)) & (time < time(nearest(time, 1.48)))); 
    indx_bsl = (time >= time(nearest(time,-0.5))) & (time < time(nearest(time, 0 ))); 
    
    % Easy condition
    tfa_easy = data.tfa.powspctrm_easy_rl(:, :, indx_all);
    tfa_easy_baseline = data.tfa.powspctrm_easy_ol(:,:,indx_bsl);
    tfa_easy_baseline = repmat(nanmean(tfa_easy_baseline, 3), [1, 1, nnz(indx_all)]);

    tfa_easy_norm = squeeze( tfa_easy ./ tfa_easy_baseline);
    log_tfa_easy = mean(log(tfa_easy_norm), 1);
                        
    % Hard condition
    tfa_hard = data.tfa.powspctrm_hard_rl(:, :, indx_all);
    tfa_hard_baseline = data.tfa.powspctrm_hard_ol(:,:,indx_bsl);
    tfa_hard_baseline = repmat(nanmean(tfa_hard_baseline, 3), [1, 1, nnz(indx_all)]);
    
    tfa_hard_norm = squeeze( tfa_hard ./ tfa_hard_baseline);    
    log_tfa_hard = mean(log(tfa_hard_norm), 1);
    
    % Prepare structures for fieldtrip tests
    ft_easy{k}.label     = {'ch1'};
    ft_easy{k}.dimord    = 'chan_freq';
    ft_easy{k}.freq      = time(indx_all);
    ft_easy{k}.powspctrm = log_tfa_easy;
    ft_easy{k} = ft_datatype_freq(ft_easy{k});  
       
    ft_hard{k}.label     = {'ch1'};
    ft_hard{k}.dimord    = 'chan_freq';
    ft_hard{k}.freq      = time(indx_all);
    ft_hard{k}.powspctrm = log_tfa_hard;
    ft_hard{k} = ft_datatype_freq(ft_hard{k});
                                 
    clear aux

end

% Statistics
cfg = [];
cfg.method      = 'montecarlo';
cfg.statistic   = 'depsamplesT';
cfg.parameter   = 'powspctrm';
cfg.numrandomization = 10000;
cfg.alpha        = 0.05;
cfg.correcttail  = 'no';
cfg.tail         = 1;
cfg.correctm    = 'cluster';
cfg.clusterstatistic='maxsum';
cfg.clusteralpha = 0.05;
cfg.clusterthreshold='parametric';

cfg.uvar=1;
cfg.ivar=2;
cfg.design(1,:) = repmat(1:n_subjects,1,2);
cfg.design(2,:) = [repmat(1,1,n_subjects) repmat(2,1,n_subjects)];

cfg.neighbours = [];

stat = ft_freqstatistics(cfg, ft_hard{:}, ft_easy{:});


% Plot average power spectrum per condition
avg_easy = ft_freqgrandaverage({}, ft_easy{:});
avg_hard = ft_freqgrandaverage({}, ft_hard{:});

figure
hold on
plot(avg_hard.freq, squeeze(avg_hard.powspctrm),'color',[0.5240,0.214,0.586],'linewidth',2.5); %% hard
plot(avg_easy.freq, squeeze(avg_easy.powspctrm),'color',[0.464,0.784,0.556],'linewidth',2.5); %% easy
plot([-1.103,-1.103], ylim,'color',[0.5240,0.214,0.586]) %hard
plot([-0.771,-0.771], ylim,'color',[0.464,0.784,0.556]) %easy
plot([0,0], ylim,'--k')
plot(xlim, [0,0],'--k')
set(gca,'fontweight','bold')

legend({'HARD','EASY'})
hold on

% Plot significative clusters
if isfield(stat, 'posclusterslabelmat')
clusters = nonzeros(unique(stat.posclusterslabelmat));
for i= 1:length(clusters)
   ix = find(stat.posclusterslabelmat == clusters(i) & stat.mask);
   if not(isempty(ix))
       plot([min(stat.freq(ix)) max(stat.freq(ix))], ...
            [0.05 * diff(ylim), 0.05 * diff(ylim)],'k-')
       [min(stat.freq(ix)) max(stat.freq(ix))]
   end
end
end

if isfield(sx,'negclusterslabelmat')
clusters = nonzeros(unique(stat.negclusterslabelmat));
for i= 1:length(clusters)
   ix = find(stat.negclusterslabelmat == clusters(i) & stat.mask);
   if not(isempty(ix))
       plot([min(stat.freq(ix)) max(stat.freq(ix))],[0.05*diff(ylim),0.05*diff(ylim)],'k-')
       [min(stat.freq(ix)) max(stat.freq(ix))]
   end
end
end

plot([-1.103, -1.103], ylim, 'color', [0.5240,0.214,0.586]) %hard
plot([-0.771, -0.771], ylim, 'color', [0.464,0.784,0.556]) %easy
plot([0, 0], ylim, '--k')
plot(xlim, [0, 0], '--k')

set(gca, 'fontweight', 'bold')