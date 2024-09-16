%% Clean-up
clear; 
close all;
% warning('off', 'MATLAB:print:FigureTooLargeForPage')
% Specify parameters
Processing2.directory = 'D:\JG\Project Gi-DREADD\Auswertung\P11\';
Processing2.timestamp = datestr(now, 'yymmdd_HHMMSS');
Processing2.age = 'P11'; % specify age group
% Get file list
Processing2.FileList = dir([Processing2.directory, '*P11', '*.mat']);
N_groups = size(Processing2.FileList, 1);

disp([num2str(N_groups), ' groups were found ...'])
% animal-wise computations
cond_all = [];
for x = 1:N_groups
    groupname = Processing2.FileList(x).name(1:end-18);
    Processing2.(groupname) = load([Processing2.directory, Processing2.FileList(x).name]); %load Processing1 variables for that animal    
    grouplabels{x} = groupname;
    group_color = Processing2.(groupname).cond_color;
    cond_all = [cond_all Processing2.(groupname).cond2];
end

%% Plots for all Groups
% occurence frequency SPWs
figure
hold on
for x = 1:N_groups
    S{x} = subplot(1,3,x);
    groupname = grouplabels{x};
    customplot_paired(Processing2.(groupname).Analysis.SPW_timewindow, Processing2.(groupname).cond2);
    ax = gca;
    set(ax.Children, 'Color', group_color(x,:))
    set(ax.Children(end, 1),'MarkerFaceColor', group_color(x,:))
    ylabel('SPW occurence (1/min)');
    title(grouplabels{x}, 'Interpreter', 'none');
end
hold off
linkaxes([S{:}],'y');
export_fig(([Processing2.directory, Processing2.age, '_', Processing2.timestamp]), '-pdf', '-painters', '-silent', '-transparent', '-nocrop', '-append'); % '-noinvert'

% progression of occurence frequency SPWs
figure
hold on
for x = 1:N_groups
    groupname = grouplabels{x};
    H = shadedErrorBar([1:1:length(Processing2.(groupname).cond_SPW)], Processing2.(groupname).Analysis.SPWprogression_mean, Processing2.(groupname).Analysis.SPWprogression_sem);
    set(H.mainLine, 'Color', group_color(x,:))
    set(H.mainLine, 'LineWidth', 2)
    set(H.patch, 'FaceColor', group_color(x,:))
    set(H.edge, 'Color', group_color(x,:))
end
hold off
legend(grouplabels, 'Interpreter', 'none')
ylabel('SPW occurence (1/min)');
title('Occurence of SPWs');
export_fig(([Processing2.directory, Processing2.age, '_', Processing2.timestamp]), '-pdf', '-painters', '-silent', '-transparent', '-nocrop', '-append'); % '-noinvert'

% occurence frequency of transient power increases
% for 4-15 Hz, which means row 4 in LFPbursts_timewindow{4,1}
figure
hold on
for x = 1:N_groups
    S{x} = subplot(1,3,x);
    groupname = grouplabels{x}; 
    customplot_paired(Processing2.(groupname).Analysis.LFPbursts_timewindow{4,1}, Processing2.(groupname).cond2);
    ax = gca;
    set(ax.Children, 'Color', group_color(x,:))
    set(ax.Children(end, 1),'MarkerFaceColor', group_color(x,:))
    ylabel('LFP bursts (1/min)');
    title(grouplabels{x}, 'Interpreter', 'none');
end
hold off
linkaxes([S{:}],'y');
ylim([0 7])
export_fig(([Processing2.directory, Processing2.age, '_', Processing2.timestamp]), '-pdf', '-painters', '-silent', '-transparent', '-nocrop', '-append'); % '-noinvert'

%% normalized Bandpower for 8-40 Hz
figure
hold on
for x = 1:N_groups
    S{x} = subplot(1,3,x);
    groupname = grouplabels{x};
    % only plotting Frequency bands 7 = 8-40 Hz
    customplot_paired(Processing2.(groupname).Analysis.bp_normalized{7, 1}, Processing2.(groupname).cond2);
    ax = gca;
    set(ax.Children, 'Color', group_color(x,:))
    set(ax.Children(end, 1),'MarkerFaceColor', group_color(x,:))
    set(ax, 'Yscale', 'log')
    ylabel('normalized Bandpower 8-40 Hz');
    title(grouplabels{x}, 'Interpreter', 'none');
    yline(1, '--r')
end
hold off
linkaxes([S{:}],'y');
export_fig(([Processing2.directory, Processing2.age, '_', Processing2.timestamp]), '-pdf', '-painters', '-silent', '-transparent', '-nocrop', '-append'); % '-noinvert'
%%
% progression of raw bandpower of subsections
figure
hold on
for x = 1:N_groups
    groupname = grouplabels{x};
    H = shadedErrorBar([1:1:length(Processing2.(groupname).condsubsec_bp)], ...
        Processing2.(groupname).Analysis.subsec_bp.all_mean, Processing2.(groupname).Analysis.subsec_bp.all_sem);
    set(H.mainLine, 'Color', group_color(x,:))
    set(H.mainLine, 'LineWidth', 2)
    set(H.patch, 'FaceColor', group_color(x,:))
    set(H.edge, 'Color', group_color(x,:))
end
hold off
legend(grouplabels, 'Interpreter', 'none')

ylabel('raw Bandpower 8-40 Hz');
title('Bandpower for subsections');
export_fig(([Processing2.directory, Processing2.age, '_', Processing2.timestamp]), '-pdf', '-painters', '-silent', '-transparent', '-nocrop', '-append'); % '-noinvert'

%% progression of normalized bandpower of subsections
figure
hold on
for x = 1:N_groups
    groupname = grouplabels{x};
    H = shadedErrorBar([1:1:length(Processing2.(groupname).condsubsec_bp)], ...
        Processing2.(groupname).Analysis.subsec_bp.all_normalized_mean, Processing2.(groupname).Analysis.subsec_bp.all_normalized_sem);
    set(H.mainLine, 'Color', group_color(x,:))
    set(H.mainLine, 'LineWidth', 2)
    set(H.patch, 'FaceColor', group_color(x,:))
    set(H.edge, 'Color', group_color(x,:))
end
legend(grouplabels, 'Interpreter', 'none')
yline(1, '--')
hold off
set(gca, 'Yscale', 'log')
ylim([0.7 1000])
ylabel('normalized Bandpower 8-40 Hz');
title('normalized Bandpower for subsections');
export_fig(([Processing2.directory, Processing2.age, '_', Processing2.timestamp]), '-pdf', '-painters', '-silent', '-transparent', '-nocrop', '-append'); % '-noinvert'

%% saving
clear H ax x S groupname
save([Processing2.directory, Processing2.age, '_', Processing2.timestamp, '.mat'], '-v7.3', '-nocompression');
disp(['File was saved to ', ...
    Processing2.directory, Processing2.age, '_', Processing2.timestamp, '.mat'])