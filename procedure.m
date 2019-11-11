% Two sample ttest2 between CB and HC at both bl and fu
for i = 1:26
    [h_tt_bl(i),p_tt_bl(i)] = ttest2(data(i,1,1:17),data(i,1,18:34));
    [h_tt_fu(i),p_tt_fu(i)] = ttest2(data(i,2,1:17),data(i,2,18:34));
end

% Two sample ttest2 between BL and FU for both HC and CB
for i = 1:26
    [h_tt_hc(i),p_tt_hc(i)] = ttest2(data(i,1,1:17),data(i,2,1:17));
    [h_tt_cb(i),p_tt_cb(i)] = ttest2(data(i,1,18:34),data(i,2,18:34));
end

% Pearson correlations
for i = 1:26
    [rho_corr_bl(i),p_corr_bl(i)] = corr(squeeze(data(i,1,18:34)),participants(18:34,1));
    [rho_corr_fu(i),p_corr_fu(i)] = corr(squeeze(data(i,2,18:34)),participants(18:34,2));
end

% Pearson correlations
for i = 1:26
    [rho_corr_usg_bl(i),p_corr_usg_bl(i)] = corr(squeeze(data(i,1,18:34)),cum_grams_before(18:34));
    [rho_corr_usg_fu(i),p_corr_usg_fu(i)] = corr(squeeze(data(i,2,18:34)),cum_grams_after(18:34));
end

% ANOVA
p_anova = zeros(26,3);
tbl_anova = cell(26,6,6);
for i = 1:26
    anova_mat = zeros(34,2);
    anova_mat(1:17,1) = data(i,1,1:17); % T1 HC
    anova_mat(18:34,1) = data(i,2,1:17); % T2 HC
    anova_mat(1:17,2) = data(i,1,18:34); % T1 CB
    anova_mat(18:34,2) = data(i,2,18:34); % T2 CB
    
    [p_anova_tmp,tbl_anova_tmp] = anova2(anova_mat,17);
    p_anova(i,:) = p_anova_tmp;
    tbl_anova{i} = tbl_anova_tmp;
    %size(p_anova)
end

% Scatter plot of volumes
reg = 7;
for i = 0:1
    y1 = squeeze(data(reg,1,(20*i+1):(20*i+20)))';
    y2 = squeeze(data(reg,2,(20*i+1):(20*i+20)))';
    d = y2 - y1;
    x = linspace(1,20,20);
    col = linspace(1,1,length(x));
    figure;
    scatter(x,y1,[],col,'filled');
    hold on;
    col = linspace(2,2,length(x));
    scatter(x,y2,[],col,'filled');
    figure;
    col = linspace(3,3,length(x));
    scatter(x,d,[],col,'filled');
end

% Arrow plot
names = {'BL'; 'FU'};
reg = 1;
y1 = squeeze(data(reg,1,21:40))';
y2 = squeeze(data(reg,2,21:40))';
y = horzcat(y1,y2);
x1 = ones(1,20);
x2 = ones(1,20)+1;
x = horzcat(x1,x2);
figure;
title('Region 1 CB');
scatter(x,y,'filled');
set(gca,'xtick',[1:2],'xticklabel',names)
for i = 1:20
    obj2 = Annotate(gca, 'line', [x1(i),x2(i)], [y1(i),y2(i)]); 
end

% Difference plot
reg = 7;
for i = 0:1
    y1 = squeeze(data(reg,1,(20*i+1):(20*i+20)))';
    y2 = squeeze(data(reg,2,(20*i+1):(20*i+20)))';
    d = y2 - y1;
    x = linspace(1,20,20);
    col = linspace(i+1,i+1,length(x));
    scatter(x,d,[],col,'filled');
    hold on;
end

% Percentage changes
for i = 1:26
    avg_vol_bl_cb(i) = mean(mean(data(i,1,18:34)));
    avg_vol_bl_hc(i) = mean(mean(data(i,1,1:17)));
    perc_bl(i) = ((avg_vol_bl_hc(i)-avg_vol_bl_cb(i))/avg_vol_bl_hc(i))*100;

    avg_vol_fu_cb(i) = mean(mean(data(i,2,18:34)));
    avg_vol_fu_hc(i) = mean(mean(data(i,2,1:17)));
    perc_fu(i) = ((avg_vol_fu_hc(i)-avg_vol_fu_cb(i))/avg_vol_fu_hc(i))*100;
end

% Cumulative usage in grams
cum_grams = zeros(34,1);
cum_grams_before = zeros(34,1);
cum_grams_between = zeros(34,1);
for i = 18:34
    cum_grams(i) = 0.5*3*52*(participants(i,1)+participants(i,2));
    cum_grams(i) = cum_grams(i) + 52*(age_at_bl(i)-age_freq_use(i))*participants(i,1);
    %
    cum_grams_before(i) = 52*(age_at_bl(i)-age_freq_use(i))*participants(i,1);
    cum_grams_between(i) = 0.5*3*52*(participants(i,1)+participants(i,2));
    cum_grams_after(i) = cum_grams_before(i) + cum_grams_between(i);
end

% MANCOVAN
Y = zeros(80,26);
groups = zeros(80,2); % bl/fu, cb/hc
covariates = zeros(80,4); % age, gender, cum_usage, cudit_bl

for i = 1:2
    for j = 1:40
        Y(40*(i-1)+j,:) = data(:,i,j)';
    end
end

for i = 1:40
    groups(i,1) = 0; % bl first 40
    groups(40+i,1) = 1; % fu next 40
end

for i = 1:20
    groups(20+i,2) = 1; % 21:40 are cb
    groups(60+i,2) = 1; % 61:80 also cb
    % 1:20, 41:60 are hc
end
groups = int16(groups);

for i = 1:40 % Ages, genders, cum_usage, cudit
    covariates(i,1) = age_at_bl(i); % bl
    covariates(40+i,1) = age_at_bl(i) + 3; % fu
    
    covariates(i,2) = genders(i);
    covariates(40+i,2) = genders(i);
    
    covariates(i,3) = cum_grams_before(i);
    covariates(40+i,3) = cum_grams_after(i);
    
    covariates(i,4) = participants(i,1);
    covariates(40+i,4) = participants(i,2);
end

[ T, p, FANCOVAN, pANCOVAN, stats ] = mancovan(Y, groups, covariates);

% Entire data matrix for use in spss

%spss_mat = zeros(80, 26+2+4); % 26 vols, 2 group vars, 4 covariates
%spss_mat(:,1:26) = Y(:,:);
%spss_mat(:,27:28) = groups(:,:);
%spss_mat(:,29:32) = covariates(:,:);
%xlswrite('dataset.xlsx',spss_mat);

spss_mat = zeros(160, 13+1+2+3); % 13 vols, 1 hemispheres, 2 group vars, age gender cum_grams as covariates

spss_mat(1:80,1:13) = Y(:,1:13); 
spss_mat(81:160,1:13) = Y(:,14:26);

spss_mat(1:80,14) = zeros(80,1);
spss_mat(81:160,14) = ones(80,1);

spss_mat(1:80,15:16) = groups(:,:);
spss_mat(81:160,15:16) = groups(:,:);

spss_mat(1:80,17:19) = covariates(:,1:3);
spss_mat(81:160,17:19) = covariates(:,1:3);

xlswrite('new_dataset.xlsx',spss_mat);

newmat = zeros(40, 52+1+3); % t1_l, t1_r, t2_l, t2_r
newmat(:,1:13) = spss_mat(1:40,1:13);
newmat(:,14:26) = spss_mat(81:120,1:13);
newmat(:,27:39) = spss_mat(41:80,1:13);
newmat(:,40:52) = spss_mat(121:160,1:13);

newmat(1:20,53) = zeros(20,1);
newmat(21:40,53) = ones(20,1);

newmat(:,54) = spss_mat(1:40,17);
newmat(:,55) = spss_mat(1:40,18);
newmat(1:20,56) = zeros(20,1);
newmat(21:40,56) = covariates(21:40,3);

xlswrite('newmat.xlsx',newmat);

median_bl = median(newmat(21:40,56));
newmat_median = newmat(21:40,:);
for i = 21:40
    if newmat(i,56) > median_bl
        newmat_median(i-20,53) = 1;
    else
        newmat_median(i-20,53) = 0;
    end
end

xlswrite('newmat_median.xlsx',newmat_median);

median_bl = median(cum_grams_before(21:40));
spss_median = zeros(160,13+1+2+3+1); % Additional median split for CBs
spss_median(:,1:19) = spss_mat(:,:);
for i = 1:160
    if spss_median(i,16) == 1
        if spss_median(i,15) == 0
            if spss_median(i,19) > median_bl
                spss_median(i,20) = 1;
            end
        else
            if spss_median(i-40,19) > median_bl
                spss_median(i,20) = 1;
            end
        end
    end
end

xlswrite('new_dataset_median.xlsx',spss_median);

%Correlations between lifetime usage stats and volumes
tmp_data = squeeze(data(:,1,21:40));
[rho_bl,pval_bl] = corr(cum_grams_bl(21:40),tmp_data');
tmp_data = squeeze(data(:,2,21:40));
[rho_fu,pval_fu] = corr(cum_grams_fu(21:40),tmp_data');

% participantwise percentage vol reductions and box plots
perc_red_sep = zeros(34,26);
for i = 1:34
    for j = 1:26
        perc_red_sep(i,j) = ((newdata(j,1,i) - newdata(j,2,i))/newdata(j,1,i))*100.0;
    end
end

y = zeros(34,26);
for i = 1:13
    y(:,2*i-1) = perc_red_sep(:,i);
    y(:,2*i) = perc_red_sep(:,13+i);
end

y_hc = y(1:17,:);
y_cb = y(18:34,:);

region_names = ["Hippocampal tail","Subiculum","CA1","Hippocampal fissure","Presubiculum", "Parasubiculum","Molecular layer","GC-ML DG","CA3","CA4","FIMBRIA","HATA","Whole Hippocampus"];

for i = 1:13
    region = i;
    figure
    %subplot(2,1,1);
    %tmp = boxplot(y_hc(:,2*region-1:2*region));
    %title(strcat('Region ',int2str(region),' HC'))
    %subplot(2,1,2);
    %tmp = boxplot(y_cb(:,2*region-1:2*region));
    %title(strcat('Region ',int2str(region),' CB'))
    
    boxplot([y_hc(:,2*region-1:2*region) y_cb(:,2*region-1:2*region)]);
    title(region_names(region));
    xticklabels({"HC_LH","HC_RH","CB_LH","CB_RH"});
    fig = gcf;
    tmp_name = strcat('volred_',int2str(i),'.jpg');
    saveas(fig, tmp_name);
end

% boxplot of hemisphere volume differences (left - right)
y = zeros(20,4,13); % 4 = t1hc, t1cb, t2hc, t2cb
for i = 1:13
    y(:,1,i) = squeeze(data(i,1,1:20) - data(i+13,1,1:20))';
    y(:,2,i) = squeeze(data(i,1,21:40) - data(i+13,1,21:40))';
    y(:,3,i) = squeeze(data(i,2,1:20) - data(i+13,2,1:20))';
    y(:,4,i) = squeeze(data(i,2,21:40) - data(i+13,2,21:40))';
end

for region = 1:13
    figure
    boxplot(y(:,:,region))
    title(region_names(region));
    xticklabels({"BL_HC","BL_CB","FU_HC","FU_CB"});
    fig = gcf;
    tmp_name = strcat('hemidiff_',int2str(region),'.jpg');
    saveas(fig,tmp_name);
end

% Median Split
lower_list = [];
higher_list = [];
median_bl = median(cum_grams_before(21:40));
for i = 21:40
    if cum_grams_before(i) < median_bl
        %append(lower_list,i-20);
        lower_list = [lower_list i-20];
    else
        %append(higher_list,i-20);
        higher_list = [higher_list i-20];
    end
end

data_median = zeros(26,2,20);

data_median(:,:,1:10) = data(:,:,lower_list+20)
data_median(:,:,11:20) = data(:,:,higher_list+20)


% SPSS dataset

spss_mat = zeros(34, 1+13*4+1+2);
for i = 1:34
    % sub id
    spss_mat(i, 1) = i;
    
    % volumes
    for v = 1:13
        spss_mat(i,1+4*v-3) = newdata(v, 1, i); % t=0, h=0
        spss_mat(i,1+4*v-2) = newdata(v+13, 1, i); % t=0, h=1
        spss_mat(i,1+4*v-1) = newdata(v, 2, i);% t=1, h=0
        spss_mat(i,1+4*v) = newdata(v+13, 2, i);% t=1, h=1
    end
    
    % group
    spss_mat(i, 58) = floor(i/18);
    
    % age_at_bl
    spss_mat(i, 59) = age_at_bl(i);
    
    % gender
    spss_mat(i, 19) = genders(i);
end

csvwrite('spss_mat_new.csv', spss_mat);

% new box plots
perc_red_sep = zeros(34,26);
for i = 1:34
    for j = 1:26
        perc_red_sep(i,j) = ((newdata(j,1,i) - newdata(j,2,i))/newdata(j,1,i))*100.0;
    end
end

y = zeros(13, 17, 4);
for v = 1:13
    y(v,:,1) = perc_red_sep(1:17,v);% HC_LH
    y(v,:,2) = perc_red_sep(1:17,v+13);% HC_RH
    y(v,:,3) = perc_red_sep(18:34,v);% CB_LH
    y(v,:,4) = perc_red_sep(18:34,v+13);% CB_RH
end

region_names = ["Hippocampal tail","Subiculum","CA1","Hippocampal fissure","Presubiculum", "Parasubiculum","Molecular layer","GC-ML DG","CA3","CA4","FIMBRIA","HATA","Whole Hippocampus"];

%for v = 1:13
%    temp = squeeze(y(v,:,:));
%    figure
%    boxplot(temp)
%    title(region_names(v));
%    xticklabels({"HC_LH","HC_RH","CB_LH","CB_RH"});
%    fig = gcf;
%    tmp_name = strcat('volred_new_',int2str(v),'.jpg');
%    saveas(fig,tmp_name);
%end

%figure
%for v = 1:13
%    temp = squeeze(y(v,:,:));
%    subplot(13,1,v);
%    boxplot(temp)
%    title(region_names(v));
%    xticklabels({"HC_LH","HC_RH","CB_LH","CB_RH"});
%end
%fig = gcf;
%tmp_name = 'volred_new_all.jpg';
%saveas(fig,tmp_name);

for v = 1:13
    temp = squeeze(y(v,:,:));
    figure
    x = horzcat(ones(1,17)-0.5,ones(1,17)+0.5,ones(1,17)+1.5,ones(1,17)+2.5);
    y1 = horzcat(y(v,:,1),y(v,:,2),y(v,:,3),y(v,:,4));
    c = horzcat(ones(1,17)+1,ones(1,17)+2,ones(1,17)+3,ones(1,17)+4);
    %boxplot(temp)
    scatter(x,y1,[],c);
    title(region_names(v));
    xticks([0 0.5 1.5 2.5 3.5]);
    xticklabels({".","HC-LH","HC-RH","CB-LH","CB-RH"});
    xlim([0 4]);
    fig = gcf;
    tmp_name = strcat('volred_scatter_new_',int2str(v),'.jpg');
    saveas(fig,tmp_name);
end