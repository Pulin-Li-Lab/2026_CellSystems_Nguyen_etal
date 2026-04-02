%% plot Fig 2C + D

% sum_all = sum(freq_all,2);
% sum_infect = sum(freq_infect,2);
% sum_ifn = sum(freq_ifn,2);
% % find mean of cell type frequency per region 
% mean_all = mean(freq_all./sum_all, 1)';
% mean_infect = mean(freq_infect./sum_infect,1, 'omitnan')';
% mean_ifn = mean(freq_ifn./sum_ifn,1, 'omitnan')';
% std_infect = std(freq_infect./sum_infect, [], 1, 'omitnan')'; 
% std_ifn = std(freq_ifn./sum_ifn, [], 1, 'omitnan')'; 
% std_all = std(freq_all./sum_all, [], 1)'; 
%%

custom_hex = ["#2CA02C", "#98DF8A", "#7F7F7F", "#8C564B", "#FFBB78", "#FF7F0E",...
    "#843C39", "#E377C2", "#AEC7E8", "#1F77B4", "#17BECF", "#BCBD22", "#C7C7C7",...
    "#F7B6D2", "#393B79", "#DBDB8D", "#9EDAE5", "#9467BD", "#FF9896", "#D62728", ...
    "#637939",  "#B5CF6B", "#C49C94"];
figure
errorbar(mean_all(:),mean_infect(:),std_infect/sqrt(10),'k')
hold on
for k = 1:length(mean_infect)
    plot(mean_all(k,:),  mean_infect(k,:),'.', 'Color',custom_hex(k), 'MarkerSize', 45)
    hold on
end
%%

figure 
errorbar(mean_infect(:),mean_ifn(:),std_ifn/sqrt(10),'k')
hold on
for k = 1:length(mean_infect)
    plot(mean_infect(k,:),  mean_ifn(k,:),'.', 'Color',custom_hex(k), 'MarkerSize', 45)
    hold on
end
ylim([0 0.12])

%% Plot Fig 2F
% separate airway and alveolar foci based on the frequency of airway epi 
% first two columns are airway

keep_ind = find(comp_foci(:,1) + comp_foci(:,2) + comp_foci(:,3)+comp_foci(:,4) ~=0)
comp = comp_foci(keep_ind,:);
alv_ind = find(comp(:,1) + comp(:,2) == 0);
epi_ind = find(comp(:,1) + comp(:,2) ~=0); 

frac = frac_l;
%plot all data usinyg scatter plot
fracp = frac(keep_ind,:);
alv_foci = fracp(alv_ind,:);
figure; scatter([1:23], fracp(alv_ind,:), 'ko')
hold on
scatter([1:23], fracp(epi_ind,:), 'kx')

% plot error bars
mean_frac = nanmean(fracp);
stderr = nanstd(fracp)/sqrt(87);
hold on; errorbar( mean_frac, stderr)

%% intra-celltype statistics on fraction of IFN+ 
epi_stat = kruskalwallis(frac(:, 1:4));
fib_stat = kruskalwallis(frac(:, 5:9));
endo_stat = kruskalwallis(frac(:, 10:13));
immune_stat = kruskalwallis(frac(:, 14:23));

%% inter-celltype statistics on fraction of IFN+ 
epi = frac(:,1:4);
fib = frac(:, 5:9);
endo = frac(:, 10:13);
immune = frac(:, 14:23);

allcell = cat(1, epi(:), fib(:), endo(:), immune(:));
g1 = repmat({'Epi'},length(epi(:)),1);
g2 = repmat({'Fib'},length(fib(:)),1);
g3 = repmat({'endo'},length(endo(:)),1);
g4 = repmat({'immune'},length(immune(:)),1);
g = [g1;g2;g3;g4];
[p,tbl,stats] = kruskalwallis(allcell, g)


%% inter-celltype statistics on fraction of IFN+ 
epi = frac(:,1:4);
fib = frac(:, 5:9);
endo = frac(:, 10:13);
immune = frac(:, 14:23);

alpha = 0.05;
numTests = 7;
p_values = zeros(numTests, 1);
[~, p_values(1)] = ttest2(epi(:), fib(:), 'Tail', 'both'); 
[~, p_values(2)] = ttest2(fib(:), endo(:), 'Tail', 'both');
[~, p_values(3)] = ttest2(endo(:), immune(:), 'Tail', 'both');
[~, p_values(4)] = ttest2(fib(:), immune(:), 'Tail', 'both');
[~, p_values(5)] = ttest2(epi(:), immune(:), 'Tail', 'both');
[~, p_values(6)] = ttest2(epi(:), endo(:), 'Tail', 'both');
[~, p_values(7)] = ttest2(fib(:), endo(:), 'Tail', 'both');

% Apply Bonferroni Correction ---
bonferroni_alpha = alpha / numTests;
fprintf('Original Alpha: %.4f, Bonferroni Alpha: %.4f\n', alpha, bonferroni_alpha);

% Determine Significant Results ---
significant_results = p_values < bonferroni_alpha;%%
%%
foci = cat(1, frac_airway, frac_alveoli);
g1 = repmat({'a'},length(frac_airway),1);
g3 = repmat({'b'},length(frac_alveoli),1);
g = [g1;g3];
g = categorical(g);
figure; boxplot(foci,g)
[p,tbl,stats] = ttest2(foci, g)

%%
foci = cat(1, star_epi, star_stromal);
g1 = repmat({'a'},length(star_epi),1);
g3 = repmat({'b'},length(star_stromal),1);
g = [g1;g3];
g = categorical(g);
figure; boxplot(foci,g)
[p,tbl,stats] = ttest2(foci, g)


