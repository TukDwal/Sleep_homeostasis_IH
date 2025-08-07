
%% MAIN
% Code supporting the findings of the manuscript "Dysregulation of slow
% wave activity in Idiopathic Hypersomnia", Tugdual Adam, Lucie Barateau,
% Yves Dauvilliers
% I am not currently able to provide the raw data which would allow to
% reproduce the figures from the article, however readers may review the
% code and apply it to their data as they wish.
% For any inquiry, contact me at : tugdual.adam@gmail.com

%% Prepare environment
clear

rng(1)

% Associated functions
% fit_S_dual_stages.m
% generateHypno.m
% finddatagroups.m
% mtspecgramc.m
% swa_ode.m
path_functions = ''; %path to associated functions
addpath(path_functions)
% Set colors for plot
c_hi1 = [17 11 92]/255;
c_hi2 = [0.4667    0.6745    0.1882];
c_hi_nlst = [0.9294    0.6941    0.1255];
c_tem = [198 8 86]/255;



%% LOAD PIPELINE
% Although I cannot presently provide the raw data, I left an example of
% the edf import pipeline, which also contains the artefact removal steps

path_edf = ''; % theoretical path to the edf folder
% Folder contains subfolders in the form Name_Surname_DDMMYY 
files = {dir(path_edf).name};

for ii = find(id_HI_nLST)'
    try
    tic
    ii
    nom = lower(br_list.NOM{ii});
    nom = replace(lower(nom),'é','e');
    nom = replace(lower(nom),'è','e');
    nom = replace(lower(nom),'â','a');
    nom = strsplit(nom,{' ','-'});
    prenom = lower(br_list.PRENOM{ii});
    prenom = replace(lower(prenom),'é','e');
    prenom = replace(lower(prenom),'è','e');
    prenom = strsplit(prenom,{' ','-'});
  
    datevisite = br_list.datevisite(ii);
    tmp1=num2str(day(datevisite));if length(tmp1)==1;tmp1=strcat('0',tmp1);end
    tmp2=num2str(month(datevisite));if length(tmp2)==1;tmp2=strcat('0',tmp2);end
    tmp3=num2str(year(datevisite)); tmp3=tmp3(3:4);
    datepattern = strcat(tmp1,tmp2,tmp3);
    
    idx = find(contains(lower(files),nom) & contains(lower(files),prenom));
    if length(idx)>1
       idtmp = find(contains(lower(files(idx)),datepattern)); 
       idx = idx(idtmp);
    end
        
    filename = files{idx};
    fold1 = fullfile(path_edf,filename);
    subfiles = {dir(fold1).name};
    edf = subfiles(find(contains(lower(subfiles),'edf')));
    hypno = subfiles(find(contains(lower(subfiles),'hypno')));
    
    % hypno: 1 = V
    %   2 = SP
    %   3 = N1
    %   4 = N2
    %   5 = N3
    
    hypno_tot = [];
    emg_eog_ecg = [];
    eeg_tot = [];
    datetime_tot = [];
    timehypno_tot = [];
    labels = {};

    
    for jj = 1:length(hypno)
        
        tmphypno = readtable(fullfile(fold1,hypno{jj}));
        
        [tmplabel,tmpdata] = edfread(fullfile(fold1,edf{jj}));
        
        if width(tmphypno)==2
            starttime = duration(str2double(strsplit(tmplabel.starttime,'.')));
            timehyp = starttime + duration(0,0,((1:height(tmphypno))-1)*30);
            timehyp = mod(timehyp,duration(24,0,0));
            tmphypno.absolutePosition_hh_mm_ss_ms_ = timehyp';
            
            tmphypno.DefaultStagingSet__stage__ = 6-tmphypno{:,2};
            tmphypno.DefaultStagingSet__stage__(tmphypno.DefaultStagingSet__stage__==6) = 9;
        elseif width(tmphypno)==3
            hypno_num = (strcmp(tmphypno{:,3},'V')) *1 + ...
                (strcmp(tmphypno{:,3},'SP') | strcmp(tmphypno{:,3},'REM'))*2 +...
                (strcmp(tmphypno{:,3},'N1') | strcmp(tmphypno{:,3},'S1'))*3 +...
                (strcmp(tmphypno{:,3},'N2') | strcmp(tmphypno{:,3},'S2'))*4 +...
                (strcmp(tmphypno{:,3},'N3') | strcmp(tmphypno{:,3},'S3'))*5 + ...
                (strcmp(tmphypno{:,3},'?')*9);
            tmphypno.DefaultStagingSet__stage__ = hypno_num;
        end
        
        hypno_tot = [hypno_tot;tmphypno.DefaultStagingSet__stage__];
        
        starttime = duration(str2double(strsplit(tmplabel.starttime,'.')));
        tmp = str2double(strsplit(tmplabel.startdate,'.'));
        startdate = datetime(tmp(3)+2000,tmp(2),tmp(1));
        startdatetime = startdate + starttime;
        fs = tmplabel.frequency(1);
        tmpdata = downsample(tmpdata',2);
        if fs==512
            tmpdata = downsample(tmpdata,2);
            fs=256;
        end
        datetime_tot = [datetime_tot;(startdatetime+duration(0,0,(1:height(tmpdata))/(fs/2)))'];
        
        tmp2 = tmphypno.absolutePosition_hh_mm_ss_ms_;
        id2 = find(diff(tmp2)<-duration(23,0,0))+1;
        for kk=id2'
            tmp2(kk:end) = tmp2(kk:end)+duration(24,0,0);
        end
   
        timehypno_tot = [timehypno_tot;startdate+tmp2];
        
        
        labels{jj} = tmplabel.label(find(contains(tmplabel.label,'C3','IgnoreCase',true) | ...
            contains(tmplabel.label,'C4','IgnoreCase',true) | ...
            contains(tmplabel.label,'O1','IgnoreCase',true) | ...
            contains(tmplabel.label,'O2','IgnoreCase',true) | ...
            contains(tmplabel.label,'M1','IgnoreCase',true) | ...
            contains(tmplabel.label,'A1','IgnoreCase',true) | ...
            contains(tmplabel.label,'M2','IgnoreCase',true) | ...
            contains(tmplabel.label,'A2','IgnoreCase',true)));
        
        eegchan = tmpdata(:,find(contains(tmplabel.label,'C3','IgnoreCase',true) | ...
            contains(tmplabel.label,'C4','IgnoreCase',true) | ...
            contains(tmplabel.label,'O1','IgnoreCase',true) | ...
            contains(tmplabel.label,'O2','IgnoreCase',true) | ...
            contains(tmplabel.label,'M1','IgnoreCase',true) | ...
            contains(tmplabel.label,'A1','IgnoreCase',true) | ...
            contains(tmplabel.label,'M2','IgnoreCase',true) | ...
            contains(tmplabel.label,'A2','IgnoreCase',true)));
        
        
        if jj == 1
            eeg_tot = [eeg_tot;[eegchan]];
            label_com = labels{jj};
        else
            [label_com,id1,id2] = intersect(labels{1},labels{2});
            eeg_tot = [eeg_tot(:,id1);eegchan(:,id2)];
        end
    end
    
    
    % Correct the '?' in the hypnogram
    id_correct = finddatagroups(hypno_tot,9);
    id_start = id_correct(1:2:end);
    dura = id_correct(2:2:end)-id_correct(1:2:end)+1;
    for idx=1:length(id_start)
        if dura(idx)<3
            if id_start(idx)>1
            hypno_tot(id_start(idx):id_start(idx)+dura(idx)-1) = hypno_tot(id_start(idx)-1);
            else
                hypno_tot(id_start(idx):id_start(idx)+dura(idx)-1) = hypno_tot(id_start(idx)+dura(idx));
            end
        else
            hypno_tot(id_start(idx):id_start(idx)+dura(idx)-1) = 1;
        end
    end
    
    % Filter
    %EEG
    [b,a] = butter(4,[0.5 40]/((fs/2)/2),'bandpass');
    eeg_tot = filtfilt(b,a,eeg_tot);
    
    % Create EEG derivations
    % C4-A1
    c4 = eeg_tot(:,find(contains(label_com,'c4','IgnoreCase',true)));
    a1 = eeg_tot(:,find(contains(label_com,'a1','IgnoreCase',true)));
    if isempty(a1)
        a1 = eeg_tot(:,find(contains(label_com,'m1','IgnoreCase',true)));
    end
    c4_a1 = c4-a1;
    
    % O2-A1
    o2 = eeg_tot(:,find(contains(label_com,'o2','IgnoreCase',true)));
    a1 = eeg_tot(:,find(contains(label_com,'a1','IgnoreCase',true)));
    if isempty(a1)
        a1 = eeg_tot(:,find(contains(label_com,'m1','IgnoreCase',true)));
    end
    o2_a1 = o2-a1;
    
    % C3-A2
    c3 = eeg_tot(:,find(contains(label_com,'c3','IgnoreCase',true)));
    a2 = eeg_tot(:,find(contains(label_com,'a2','IgnoreCase',true)));
    if isempty(a2)
        a2 = eeg_tot(:,find(contains(label_com,'m2','IgnoreCase',true)));
    end
    c3_a2 = c3-a2;
    

    if isempty(o2)
        o1 = eeg_tot(:,find(contains(label_com,'o1','IgnoreCase',true)));
        o2_a1 = o1-a2;
    end
    
    if isempty(o2_a1)
        o2_a1 = zeros(1,length(c4_a1));
    end
    % ARTEFACT REMOVAL
    
    % Start by creating evenly spaced data and synchronize hypnogram
   
    time_eeg = datetime_tot(1):duration(0,0,1/(fs/2)):datetime_tot(end);
    c4_a1 = interp1(datetime_tot,c4_a1,time_eeg);
    o2_a1 = interp1(datetime_tot,o2_a1,time_eeg);
    c3_a2 = interp1(datetime_tot,c3_a2,time_eeg);
    
    hypno = interp1(timehypno_tot,hypno_tot,time_eeg,'nearest');
    hypno(hypno>6) = 1;

    % Amplitude threshold on the derivative 
    % Split up SLEEP and WAKE
    
    % C4-A1
    TMP = {c4_a1;o2_a1;c3_a2};
    for jj=1:3
        % SLEEP
        data_sleep = TMP{jj};
        data_sleep(hypno==1) = NaN;
        diff_data_sleep = [0 diff(data_sleep)];
        artef = abs(diff_data_sleep)>prctile(abs(diff_data_sleep),99.9);
        artef = conv(artef,ones(1+2*256,1),'same');
        data_sleep(artef>1) = nan;
        
        % WAKE
        data_wake = TMP{jj};
        data_wake(hypno>1) = NaN;
        diff_data_wake = [0 diff(data_wake)];
        artef = abs(diff_data_wake)>prctile(abs(diff_data_wake),98);
        artef = conv(artef,ones(1+2*256,1),'same');
        data_wake(artef>1) = nan;
        
        % Reconstruct artefact free signal
        data_sleep(isnan(data_sleep)) = 0;
        data_wake(isnan(data_wake)) = 0;
        tmp = data_sleep + data_wake;
        tmp(tmp==0) = NaN;
        TMP{jj} = tmp;
    end
        
    % Spectral analysis    
    params=struct('tapers',[],'Fs',[],'fpass',[]);
    params.tapers = [5/2 4];
    params.Fs = fs/2;
    params.fpass = [0.3 40];
    movwind = [4 2];
    
    
    VEC_spec = {};
    for jj=1:3
        
        [S,t,f] = mtspecgramc(TMP{jj},movwind,params);
        S = real(S);
        idnan = isnan(S);
        [t_grid,f_grid] = meshgrid(t,f);
        
        Snan=griddata(t_grid(idnan==0),f_grid(idnan==0), S(idnan==0),...
            t_grid(idnan==1),f_grid(idnan==1),'natural');
        
        Sinterp = S;
        if sum(isnan(Sinterp(:))) == length(Sinterp(:)) % If no O2, set to 0 for this channel (otherwise error)
            Sinterp = zeros(size(Sinterp));
        else
            Sinterp(idnan==1) = Snan;
        end
        
        S_pxx = Sinterp.*f;        
        % Delta low 0.5-2Hz
        [~,id1] = min(abs(f-0.5));
        [~,id2] = min(abs(f-2));
        delta_1 = trapz(f(id1:id2),S_pxx(:,id1:id2),2);
        % Delta high 2-4Hz
        [~,id1] = min(abs(f-2));
        [~,id2] = min(abs(f-4));
        delta_2 = trapz(f(id1:id2),S_pxx(:,id1:id2),2);
        % Theta 4-7Hz
        [~,id3] = min(abs(f-7));
        theta_1 = trapz(f(id2:id3),S_pxx(:,id2:id3),2); %theta_1(pow_tag2) = nan;
        % Alpha 8-12Hz
        [~,id4] = min(abs(f-8));
        [~,id5] = min(abs(f-12));
        alpha_1 = trapz(f(id4:id5),S_pxx(:,id4:id5),2); %alpha_1(pow_tag2) = nan;
        % Sigma 12-16Hz
        [~,id6] = min(abs(f-16));
        sigma_1 = trapz(f(id5:id6),S_pxx(:,id5:id6),2); %sigma_1(pow_tag2) = nan;
        % Beta 15-25Hz
        [~,id7] = min(abs(f-16));
        [~,id8] = min(abs(f-25));
        beta_1 = trapz(f(id7:id8),S_pxx(:,id7:id8),2); %beta_1(pow_tag2) = nan;
        % Low Gamma 25-40Hz
        [~,id9] = min(abs(f-40));
        gamma_1 = trapz(f(id8:id9),S_pxx(:,id8:id9),2); %gamma_1(pow_tag2) = nan;
        
        VEC_spec{jj} = [delta_1,delta_2,theta_1,alpha_1,sigma_1,beta_1,gamma_1];
    end

    t_spectr = duration(0,0,t);
    idx = find(diff(timehypno_tot)<-duration(22,0,0))+1;
    for jj=idx'
        timehypno_tot(jj:end) = timehypno_tot(jj:end)+duration(24,0,0);
    end
    [~,idunique] = unique(timehypno_tot);
    hypno_spectr = interp1(timehypno_tot(idunique)-timehypno_tot(1),hypno_tot(idunique),t_spectr,'nearest');
    
    tab_spec = table(t_spectr',hypno_spectr',...
        VEC_spec{1}(:,1),VEC_spec{1}(:,2),VEC_spec{1}(:,3),VEC_spec{1}(:,4),VEC_spec{1}(:,5), VEC_spec{1}(:,6), VEC_spec{1}(:,7),...
        VEC_spec{2}(:,1),VEC_spec{2}(:,2),VEC_spec{2}(:,3),VEC_spec{2}(:,4),VEC_spec{2}(:,5), VEC_spec{2}(:,6), VEC_spec{2}(:,7),...
        VEC_spec{3}(:,1),VEC_spec{3}(:,2),VEC_spec{3}(:,3),VEC_spec{3}(:,4),VEC_spec{3}(:,5), VEC_spec{3}(:,6), VEC_spec{3}(:,7),...
        'VariableNames',...
        {'Time','Hypno','[0.5-2]hz C4-A1','[2-4]hz C4-A1','[4-7]hz C4-A1','[8-12]hz C4-A1','[12-16]hz C4-A1','[16-25]hz C4-A1','[25-40]hz C4-A1',...
        '[0.5-2]hz O2-A1','[2-4]hz O2-A1','[4-7]hz O2-A1','[8-12]hz O2-A1','[12-16]hz O2-A1','[16-25]hz O2-A1','[25-40]hz O2-A1',...
        '[0.5-2]hz C3-A2','[2-4]hz C3-A2','[4-7]hz C3-A2','[8-12]hz C3-A2','[12-16]hz C3-A2','[16-25]hz C3-A2','[25-40]hz C3-A2'});
    
    writetable(tab_spec,fullfile(fold1,'puissances_spectrales.xlsx'));
    spectral_tables{ii} = tab_spec; 
    end
end


%% Load processed data - to load from saved pre-processed data

spectral_tables = {};
path 
for ii=find(id_HI_LST + id_HI_nLST + id_tem)'
    ii
    nom = lower(br_list.NOM{ii});
    nom = replace(lower(nom),'é','e');
    nom = replace(lower(nom),'è','e');
    nom = replace(lower(nom),'â','a');
    nom = strsplit(nom,{' ','-'});
    prenom = lower(br_list.PRENOM{ii});
    prenom = replace(lower(prenom),'é','e');
    prenom = replace(lower(prenom),'è','e');
    prenom = strsplit(prenom,{' ','-'});
    datevisite = br_list.datevisite(ii);
    tmp1=num2str(day(datevisite));if length(tmp1)==1;tmp1=strcat('0',tmp1);end
    tmp2=num2str(month(datevisite));if length(tmp2)==1;tmp2=strcat('0',tmp2);end
    tmp3=num2str(year(datevisite)); tmp3=tmp3(3:4);
    datepattern = strcat(tmp1,tmp2,tmp3);
    
   
    subfolders = {dir('F:\PhD\BedRest\Analyse_spectral_BR').name};
    subfolders = subfolders(3:end);
    
    for kk = 1:length(subfolders)
        path = fullfile('F:\PhD\BedRest\Analyse_spectral_BR',subfolders{kk});
        subfiles = {dir(path).name};
        idx = find(contains(lower(subfiles),nom) & contains(lower(subfiles),prenom));
       
        if length(idx)>1
            idx = find(contains(lower(subfiles),nom) & contains(lower(subfiles),prenom) & contains(lower(subfiles),datepattern));
        end
        if isempty(idx)
            continue
        end
        subfiles2 = {dir(fullfile(path,subfiles{idx})).name};
        
        spec_tab = subfiles2(find(contains(subfiles2,'puissances_spectrales')));
        
        opts = detectImportOptions(fullfile(fullfile(path,subfiles{idx}),spec_tab{1}));
        opts.VariableTypes{1} = 'duration';
        datatab = readtable(fullfile(fullfile(path,subfiles{idx}),spec_tab{1}),opts);
        spectral_tables{ii} = datatab;
    end   
end


%% DETECT SLEEP CYCLES

% spectral_tables is a structure array containing 1 table for each
% participant
% The table contains the followingg columns:
% - time
% - hypnogram (5:NREM3, 4:NREM2, 3:NREM1, 2:REM, 1:W)
% - spectral power in [0.5-2]Hz, [2-4]Hz, [4-7.5]Hz, [8-12]Hz, [12-16]Hz, [16-24]Hz
% in C3-M2, C4-M2 and O2-M1
% datapoints are sampled every 2 seconds

for ii=1:length(spectral_tables)
    ii
    tab = spectral_tables{ii};
    if ~isempty(tab)
    hypno = tab.Hypno;
    time = tab.Time;
    
    % "Movmedian" of 5 min to smooth
    for kk=1:length(hypno)
        hypno(kk) = mode(hypno(max(1,kk-5*15):min(length(hypno),kk+5*15)));       
    end
    
    %Detection algorithm
    
    
    k_cycle = 0;
    id = 1;
    cycles_tag = zeros(1,length(hypno));
    ENDO = 0; % 0 if awake, 1 if asleep
    while id<length(hypno)-449
        % INIT = awake : search for the 1st epoch of sleep followed by 15min of stable 
        while ENDO == 0 & id<length(hypno)-449
            hyp15min = sum(hypno(id:id+449)==1); %15min =  15*30 = 450 datapoints  
         
            if hyp15min>0
                id=id+1;
                hyp15min = sum(hypno(id:id+449)==1);
            else
                ENDO = 1;
                k_cycle = k_cycle+1;
            end
        end
        
        % STEP 1 : sleep onset found : cycle starts, we now search for end
        % of cycle:
        %  15min of wake = end of cycle
        % If REM, look for end of REM = 15min of stable not-REM (W or NREM) 
        
        while ENDO == 1 & id<length(hypno)-449
            if hypno(id) == 1 %if wake
                % extract following 15min
                hyp15min = sum(hypno(id:id+449)==1);
                if hyp15min==450 % if stable wake, iterate cycle
                    ENDO = 0;
                    id = id+450;
                else % otherwise same cycle
                    cycles_tag(id) = k_cycle;
                    id = id+1;
                end
            elseif hypno(id) == 2 %if REM
                % extract following 15 min
  
                hyp15min = sum(hypno(id+1:id+450)==2);
                if hyp15min==0 %if no REM, iterate cycle
                    cycles_tag(id) = k_cycle;
                    k_cycle = k_cycle+1;
                    id = id+1;
                else     % otherwise same cycle
                    cycles_tag(id) = k_cycle;
                    id = id+1;
                end
                  
                % If first cycle, NREM3>NREM2 and NREm2 is stable on30min : iterate cycle 
           elseif hypno(id) == 5 & hypno(id+1) == 4 & k_cycle==1 & sum(cycles_tag==k_cycle)/30>45
                 
                hyp30min = sum(hypno(id+1:id+900)>=4);
                if hyp30min==900 %30min NREM2 = iterate cycle
                    cycles_tag(id) = k_cycle;
                    k_cycle = k_cycle+1;
                    id = id+1;
                else     % otherwise same cycle
                    cycles_tag(id) = k_cycle;
                    id = id+1;
                end
                
                
            else 
                cycles_tag(id) = k_cycle;
                id = id+1;
            end
        end
        
    end
    if sum(hypno(end-449:end)>1) >0
        cycles_tag(end-449:end) = k_cycle;
    end
    % If a cycle lasts less than 15imn, add to previous one
    for kk=1:k_cycle
        cycl = cycles_tag(cycles_tag==kk);
        if length(cycl) < 450
            cycles_tag(cycles_tag==kk) = kk-1;
        end
    end
    
    % Throw first cycle if <30min
    length1 = sum(cycles_tag==1);
    if length1<30*30
        cycles_tag(cycles_tag==1) = 0;
    end
    % Fix consecutive cycle numbers (if 1 skipped)
    uni_cycles = unique(cycles_tag);
    id_skip = find(diff(uni_cycles)==2)+1;
    for jj = id_skip
        cycles_tochange = uni_cycles(jj:end);
        for kk=cycles_tochange
            cycles_tag(cycles_tag==kk) = cycles_tag(cycles_tag==kk)-1;
            uni_cycles = unique(cycles_tag);
        end
    end
    
    cycles_tag(hypno==0) = 0;

    tab.Cycles = cycles_tag  % We add the tagged cycle to the table
    
    % Cycles = 0 means wake, or out of cycle, otherwise cycles are numbered
    % in consecutive order
    
    spectral_tables{ii} = tab;
    end
end

%% COMPUTE SWA

id_c4_m2 = [19,50,66,69,103,105,112,122,133,143,161,166,169,169,236,247,283,308,310,331,346,371,375,406,414,417,429,449,452,459,460,463,496,497,501,506,535,561,583,587,588,590,600];
% by default I am analyzing C4-M1, but for some participants, C4-M1 is
% non-exploitable, and I'm using C4-M2 instead. Those partiipants are
% listed in id_c4_m2


swa_array = {};
theta_array = {};

for ii=1:length(spectral_tables)
    ii
    tab =  spectral_tables{ii};
    if isempty(tab)
        continue
    end
    time = tab.Time;
    hypno = tab.Hypno;
 
    if ismember(ii, id_c4_m2) 
        swa = (tab{:,17}+tab{:,18});
    else
        swa = (tab{:,3}+tab{:,4});
    end
    
    swa = movmedian(swa,30,'omitnan'); % smooth to remove spurious residual artefacts
    
    % Normalize by swa over first 7 hours of sleep
    swa = swa / mean(swa(hypno<=5 & hypno>=2 & time<duration(7,0,0)),'omitnan');
    spectral_tables{ii}.swa = swa; % save in table
    
   val_swa = [];

    for jj=setdiff(unique(cycles),0)'
        if sum(cycles==jj)/30 < 45 & ~(cycles(max(1,find(cycles==jj,1,'first')-1))==0 & cycles(min(find(cycles==jj,1,'last')+1,length(cycles)))==0)
            continue
        end
        swa_cycle = swa(cycles==jj);
        hypno_cycle = hypno(cycles==jj);
        time_cycle = hours(time(cycles==jj));
        deltapow_auc = trapz(time_cycle(hypno_cycle>3),swa_cycle(hypno_cycle>3));  % AUC of delta power
        deltapow_med = median(swa(cycles==jj & hypno>3),'omitnan'); % median of SWA power in N2/N3 only
        deltapow_max = prctile(swa_cycle(hypno_cycle>3),95); % 95th percentile
        midcycle = median(time_cycle); % for time: either midcycle
        startcycle = time_cycle(1); % or start of cycle
        val_swa = [val_swa [deltapow_med;deltapow_auc;deltapow_max;midcycle;startcycle]];
    end
    val_swa(:,isnan(val_swa(1,:))) = [];
    
    [~,idsort] = sort(val_swa(5,:));
    swa_array{ii} = val_swa(:,idsort);
end


%% COMPARE DELTA C3-M1 / C4-M2

% Here we compare LIN's concordance coefficient for participants where both
% are exploitable
for ii=1:length(spectral_tables)
    ii
    tab =  spectral_tables{ii};
    if isempty(tab)
        continue
    end
    time = tab.Time;
    hypno = tab.Hypno;
    cycles = tab.Cycles;
    [~,id] = min(abs(hours(time)-[0:2:32]));
    if ismember(ii, id_c4_m2)
        continue
    else
        
        swa = (tab{:,17}+tab{:,18}); % Columns for delta C4-M2
        
    end
    

    swa = movmedian(swa,30,'omitnan'); % retire les artefacts residuels
    
    % Normalize
    swa = swa / mean(swa(hypno<=5 & hypno>=2 & time<duration(7,0,0)),'omitnan');
    spectral_tables{ii}.swa2 = swa; % save as the SWA for C4-A2 
    
end

% Correlation coff swa and swa 2 

Rtot = nan(1,619);
CCCtot = nan(619,3);
for ii=1:length(spectral_tables)
    ii
    tab =  spectral_tables{ii};
    if isempty(tab) | sum(strcmp(tab.Properties.VariableNames,'swa2'))==0 | sum(strcmp(tab.Properties.VariableNames,'swa'))==0
        continue
    end
    
    swa1 = tab.swa;
    swa2 = tab.swa2;
    hypno = tab.Hypno;
    swa1 = swa1(hypno>=3);
    swa2 = swa2(hypno>=3);
    % Remove probable artefacts
    swa1(swa1>10) = NaN;
    swa2(swa2>10) = NaN;
    R = corrcoef(swa1,swa2, 'rows','complete');
    Rtot(ii) = R(1,2);
    idnan = isnan(swa1) + isnan(swa2);
    ccc = f_CCC([swa1(~idnan),swa2(~idnan)],0.05);
    ccc = ccc{1};
    CCCtot(ii,:) = [ccc.est, ccc.confInterval];

end
% Manually removed because in those participants C3-M1 was ok but not c4-M2
Rtot([92 101 256 261 262 355 356 445 545]) = NaN;
CCCtot([92 101 256 261 262 355 356 445 545],:) = NaN;

% Distribution per diagnose
sum(~isnan(CCCtot(id_4g==1)))
sum(~isnan(CCCtot(id_4g==2)))
sum(~isnan(CCCtot(id_4g==3)))
sum(~isnan(CCCtot(id_4g==4)))

x = 0:0.01:100;
y = sum(CCCtot(find(id_4g>0),1)>x) / sum(~isnan(CCCtot(find(id_4g>0),1)));

figure
plot(x,y)
ylabel('Proportion (%)')
xlabel('CCC')
title('Proportion of observations above a CCC threshold')

%% FIGURE 2

% Proportion of sleep per hour of recording

tts_all = nan(619,32,3);

for ii=1:619
   tab = spectral_tables{ii};
   if isempty(tab)
       continue
   end
   time = tab.Time;
   hypno = tab.Hypno;
   [~,idcut] = min(abs(hours(time) - [0:1:32]));
   for jj=1:length(idcut)-1
      hypno_window = hypno(idcut(jj):idcut(jj+1)-1);
      tts_all(ii,jj,2) = sum(hypno_window>=2)*2;
   end   
end

f2a=figure
hold on
plot(1:32,mean(tts_all(id_4g==1,:,2) + tts_all(id_4g==1,:,3))/60,'color',c_hi1,'LineWidth',2)
plot(1:32,mean(tts_all(id_4g==2,:,2) + tts_all(id_4g==2,:,3))/60,'color',c_hi2,'LineWidth',2)
plot(1:32,mean(tts_all(id_4g==3,:,2) + tts_all(id_4g==3,:,3))/60,'color',c_hi_nlst,'LineWidth',2)
plot(1:32,mean(tts_all(id_4g==4,:,2) + tts_all(id_4g==4,:,3))/60,'color',c_tem,'LineWidth',2)
ylabel('Duration (min)')
title('Sleep time per hour')
xticks((1:32)-0.5)
xticklabels(string(mod(23:55,24)))
xlabel('Time (hours)')
grid on

%% FIGURE 3

% SWA energy on 2 consecutive nychthemerons
energy_d1_d2 = [];
rem_dura_d1_d2 = [];

for ii=find(id_4g>0)'
    
    tab = spectral_tables{ii};
    swa = tab.swa;
    hypno = tab.Hypno;
    time = tab.Time;
    [~,id1] = min(abs(time-duration(8,0,0))); % 8 hours in = 11pm-7am night 1
    [~,id2] = min(abs(time-duration(32,0,0)));
    swa_D1 = swa(1:id1);
    swa_D2 = swa(id1+1:end);
    
    hypno_D1 = hypno(1:id1);
    hypno_D2 = hypno(id1+1:end);
    swa_D1(hypno_D1<=2) = 0;
    swa_D2(hypno_D2<=2) = 0;
    
    swa_energy_d1 = sum(swa_D1);
    swa_energy_d2 = sum(swa_D2);
    
    rem_dura_d1 = sum(hypno_D1==2)/30;
    rem_dura_d2 = sum(hypno_D2==2)/30;
    energy_d1_d2(ii,:) = [swa_energy_d1,swa_energy_d2];
    rem_dura_d1_d2(ii,:) = [rem_dura_d1,rem_dura_d2];
    
end

f3=figure
t=tiledlayout(1,4,'TileSpacing','tight')
t1=nexttile;
hold on
energy_1 = [];
for ii=find(id_4g==1)'
    energy_d1 = energy_d1_d2(ii,1);
    x_d1 = 0.9 + rand*0.2;
    energy_d2 = energy_d1_d2(ii,2);
    x_d2 = 1.9 + rand*0.2;
    line([x_d1,x_d2],[energy_d1,energy_d2],'color',[c_hi1 0.15])
    scatter([x_d1,x_d2],[energy_d1,energy_d2],[],c_hi1,'filled','MarkerFaceAlpha',0.15)
    energy_1 = [energy_1;[energy_d1,energy_d2]];
end
line([1 2],nanmean(energy_1),'color',[c_hi1],'LineWidth',3)
s1=scatter([1 2],nanmean(energy_1),30,c_hi1,'filled')
xticks(1:2)
xticklabels({'Previous day + night 1','Day 1 + night 2'})
ylabel('SWA energy')
t2=nexttile;
hold on
energy_2 = [];
for ii=find(id_4g==2)'
    energy_d1 = energy_d1_d2(ii,1);
    x_d1 = 0.8 + rand*0.4;
    energy_d2 = energy_d1_d2(ii,2);
    x_d2 = 1.8 + rand*0.4;
    line([x_d1,x_d2],[energy_d1,energy_d2],'color',[c_hi2 0.15])
    scatter([x_d1,x_d2],[energy_d1,energy_d2],[],c_hi2,'filled','MarkerFaceAlpha',0.15)
    energy_2 = [energy_2;[energy_d1,energy_d2]];
end
line([1 2],nanmean(energy_2),'color',[c_hi2],'LineWidth',3)
s2=scatter([1 2],nanmean(energy_2),30,c_hi2,'filled')
xticks(1:2)
xticklabels({'Previous day + night 1','Day 1 + night 2'})
ylabel('SWA energy')
t3=nexttile;
hold on
energy_3 = [];
for ii=find(id_4g==3)'
    energy_d1 = energy_d1_d2(ii,1);
    x_d1 = 0.8 + rand*0.4;
    energy_d2 = energy_d1_d2(ii,2);
    x_d2 = 1.8 + rand*0.4;
    line([x_d1,x_d2],[energy_d1,energy_d2],'color',[c_hi_nlst 0.25])
    scatter([x_d1,x_d2],[energy_d1,energy_d2],[],c_hi_nlst,'filled','MarkerFaceAlpha',0.25)
    energy_3 = [energy_3;[energy_d1,energy_d2]];
end
line([1 2],nanmean(energy_3),'color',[c_hi_nlst],'LineWidth',3)
s3=scatter([1 2],nanmean(energy_3),30,c_hi_nlst,'filled')
xticks(1:2)
xticklabels({'Previous day + night 1','Day 1 + night 2'})
ylabel('SWA energy')
t4=nexttile;
hold on
energy_4 = [];
for ii=find(id_4g==4)'
    energy_d1 = energy_d1_d2(ii,1);
    x_d1 = 0.8 + rand*0.4;
    energy_d2 = energy_d1_d2(ii,2);
    x_d2 = 1.8 + rand*0.4;
    line([x_d1,x_d2],[energy_d1,energy_d2],'color',[c_tem 0.15])
    scatter([x_d1,x_d2],[energy_d1,energy_d2],[],c_tem,'filled','MarkerFaceAlpha',0.15)
    energy_4 = [energy_4;[energy_d1,energy_d2]];
end
line([1 2],nanmean(energy_4),'color',[c_tem],'LineWidth',3)
s4=scatter([1 2],nanmean(energy_4),30,c_tem,'filled')
xticks(1:2)
xticklabels({'Previous day + night 1','Day 1 + night 2'})
ylabel('SWA energy')
lgd = legend([s1 s2 s3 s4],{'IH-LST-S','IH-LST-nS','IH-nLST','HC'},'NumColumns',4)
lgd.Layout.Tile = "south"
xlim([0.75 2.25])
linkaxes([t1 t2 t3 t4],'xy')
title(t,'Comparison of SWA energy (NREM pressure) and REM duration (REM pressure) on 2 consecutive days')

% STATS
% ENERGY
% HI-LST-S
tmptbl = table(energy_1(:,1),energy_1(:,2),br_list.Age(id_4g==1),categorical(br_list.SEXE(id_4g==1)),...
    categorical([1:sum(id_4g==1)]'),'VariableNames',{'E1','E2','Age','Sex','Pat'});
mdl = fitglme(tmptbl, 'E1 ~ E2 + Age + Sex +(1|Pat)') 
p1 = mdl.Coefficients.pValue(2);
% HI-LST-nS
tmptbl = table(energy_2(:,1),energy_2(:,2),br_list.Age(id_4g==2),categorical(br_list.SEXE(id_4g==2)),...
    categorical([1:sum(id_4g==2)]'),'VariableNames',{'E1','E2','Age','Sex','Pat'});
mdl = fitglme(tmptbl, 'E1 ~ E2 + Age + Sex +(1|Pat)') 
p2 = mdl.Coefficients.pValue(2);
% HI-nLST
tmptbl = table(energy_3(:,1),energy_3(:,2),br_list.Age(id_4g==3),categorical(br_list.SEXE(id_4g==3)),...
    categorical([1:sum(id_4g==3)]'),'VariableNames',{'E1','E2','Age','Sex','Pat'});
mdl = fitglme(tmptbl, 'E1 ~ E2 + Age + Sex +(1|Pat)') 
p3 = mdl.Coefficients.pValue(2);
% TEM
tmptbl = table(energy_4(:,1),energy_4(:,2),br_list.Age(id_4g==4),categorical(br_list.SEXE(id_4g==4)),...
    categorical([1:sum(id_4g==4)]'),'VariableNames',{'E1','E2','Age','Sex','Pat'});
mdl = fitglme(tmptbl, 'E1 ~ E2 + Age + Sex +(1|Pat)') 
p4 = mdl.Coefficients.pValue(2);

p_adj = mafdr([p1,p2,p3,p4],'BHFDR',true)

%rmANOVA to compare linear evolution between groups

t_swa = table(categorical([ones(length(energy_1(:,1)),1);ones(length(energy_2(:,1)),1)*2;ones(length(energy_3(:,1)),1)*3;ones(length(energy_4(:,1)),1)*4]),...
    [energy_1(:,1);energy_2(:,1);energy_3(:,1);energy_4(:,1)],...
    [energy_1(:,2);energy_2(:,2);energy_3(:,2);energy_4(:,2)],...
    [br_list.BEDREST_TTS32(id_4g==1);br_list.BEDREST_TTS32(id_4g==2);br_list.BEDREST_TTS32(id_4g==3);br_list.BEDREST_TTS32(id_4g==4)],...
     [br_list.Age(id_4g==1);br_list.Age(id_4g==2);br_list.Age(id_4g==3);br_list.Age(id_4g==4)],...
     categorical([br_list.SEXE(id_4g==1);br_list.SEXE(id_4g==2);br_list.SEXE(id_4g==3);br_list.SEXE(id_4g==4)]),...
'VariableNames',{'group','d1','d2','tst','Age','Sex'});
rmtime = table([1 2]','VariableNames',{'rmtime'});
rm_swa = fitrm(t_swa,'d1-d2~group + Age + Sex','WithinDesign',rmtime);
ranovatbl_swa = ranova(rm_swa)

result_swa = multcompare(rm_swa,'group')


%% FIGURE 4

% SLEEP CYCLE CHARACTERISTICS

wake_up_index = nan(619,1); % first find the wake-up time = first occurrence of 1hrs of stable wake/N1
for ii=find(id_4g>0)'
    tab = spectral_tables{ii};
    hypno = tab.Hypno;
    cycles = tab.Cycles;
    id1 = find(cycles>0,1,'first');
    conv_w = conv(hypno(id1+1:end)==1 | hypno(id1+1:end)==3, ones(1,60*30),'valid');
    id2 = find(conv_w==1800,1,'first') + id1;
    if isempty(id2)     
        conv_w = conv(hypno(id1+1:end)==1 | hypno(id1+1:end)==3, ones(1,45*30),'valid'); % if never 1hour of stable wake, allow 45min instead
        id2 = find(conv_w==45*30,1,'first') + id1;
    end
    wake_up_index(ii) = id2;
end

charac_c = {};
n_cycles = [];
for ii=find(id_4g>0)'
   tab = spectral_tables{ii};
   if isempty(tab)
       continue
   end
   cycles = tab.Cycles;
   wuid = wake_up_index(ii);
   id1 = find(cycles>0,1,'first');
   cycles_n1 = cycles(id1:wuid);
   dura_c = [];
   mid_c = [];
   for j=setdiff(unique(cycles_n1),0)'
       dura_c(end+1) = sum(cycles_n1==j)/30;
       mid_c(end+1) = mean(find(cycles_n1==j));
   end
   charac_c{ii} = [dura_c;mid_c];
   n_cycles(ii,:) = [length(mid_c),(wuid-id1)/30];
end

% cycle duration across the night, align on cycles
dura_c_align = nan(619,16);
for ii=find(id_4g>0)'
   tab = charac_c{ii};
   if isempty(tab)
       continue
   end
   dura_c = charac_c{ii}(1,:);   
   mid_c = charac_c{ii}(2,:);
   dura_c_align(ii,1:length(dura_c)) = dura_c;
end

f4ab=figure
t=tiledlayout(1,2,'TileSpacing','tight')
t1=nexttile;
hold on
scatter(n_cycles(id_4g==1,2)/60,n_cycles(id_4g==1,1),[],[c_hi1],'filled','MarkerFaceAlpha',0.7)
scatter(n_cycles(id_4g==2,2)/60,n_cycles(id_4g==2,1),[],[c_hi2],'filled','MarkerFaceAlpha',0.7)
scatter(n_cycles(id_4g==3,2)/60,n_cycles(id_4g==3,1),[],[c_hi_nlst],'filled','MarkerFaceAlpha',0.7)
scatter(n_cycles(id_4g==4,2)/60,n_cycles(id_4g==4,1),[],[c_tem],'filled','MarkerFaceAlpha',0.7)
mdl1 = fitlm(n_cycles(id_4g==1,2)/60,n_cycles(id_4g==1,1))
y1 = mdl1.Coefficients.Estimate(1) + mdl1.Coefficients.Estimate(2)*(0:0.5:25);
p1=plot((0:0.5:25),y1,'color',c_hi1,'LineWidth',2)
mdl2 = fitlm(n_cycles(id_4g==2,2)/60,n_cycles(id_4g==2,1))
y2 = mdl2.Coefficients.Estimate(1) + mdl2.Coefficients.Estimate(2)*(0:0.5:25);
p2=plot((0:0.5:25),y2,'color',c_hi2,'LineWidth',2)
mdl3 = fitlm(n_cycles(id_4g==3,2)/60,n_cycles(id_4g==3,1))
y3 = mdl3.Coefficients.Estimate(1) + mdl3.Coefficients.Estimate(2)*(0:0.5:25);
p3=plot((0:0.5:25),y3,'color',c_hi_nlst,'LineWidth',2)
mdl4 = fitlm(n_cycles(id_4g==4,2)/60,n_cycles(id_4g==4,1))
y4 = mdl4.Coefficients.Estimate(1) + mdl4.Coefficients.Estimate(2)*(0:0.5:25);
p4=plot((0:0.5:25),y4,'color',c_tem,'LineWidth',2)
xlabel('Sleep duration (hours)')
ylabel('Number of sleep cycles across the first main sleep episode')
grid on
t2=nexttile;
hold on
tmp1 = dura_c_align(id_4g==1,1:13);
id1 = sum(~isnan(tmp1))<=5;
tmp1(:,id1) = NaN;
errorbar(1:13,mean(tmp1,'omitnan'),std(tmp1,'omitnan'),'-s','color',c_hi1,"MarkerFaceColor",c_hi1)
tmp2 = dura_c_align(id_4g==2,1:13);
id2 = sum(~isnan(tmp2))<=5;
tmp2(:,id2) = NaN;
errorbar((1:13)+0.1,mean(tmp2,'omitnan'),std(tmp2,'omitnan'),'-s','color',c_hi2,"MarkerFaceColor",c_hi2)
tmp3 = dura_c_align(id_4g==3,1:13);
id3 = sum(~isnan(tmp3))<=5;
tmp3(:,id3) = NaN;
errorbar((1:13)+0.2,mean(tmp3,'omitnan'),std(tmp3,'omitnan'),'-s','color',c_hi_nlst,"MarkerFaceColor",c_hi_nlst)
tmp4 = dura_c_align(id_4g==4,1:13);
id4 = sum(~isnan(tmp4))<=5;
tmp4(:,id4) = NaN;
errorbar((1:13)+0.3,mean(tmp4,'omitnan'),std(tmp4,'omitnan'),'-s','color',c_tem,"MarkerFaceColor",c_tem)
xlabel('Consecutive cycles')
ylabel('Duration (min)')
xticks(1:13)
xticklabels(string(1:13))
grid on
title(t,'Characterization of cycles across the different groups')
lgd = legend([p1 p2 p3 p4],{'IH-LST-S','IH-LST-nS','IH-nLST','HC'},'NumColumns',4)
lgd.Layout.Tile = "south"

% STATISTICS 

% ANCOVA on linear fits

[h,atab,ctab,stats] = aoctool(n_cycles(id_4g>0,2)/60,n_cycles(id_4g>0,1),id_4g(id_4g>0));

% GLMs across consecutive cycles
pval_c = [];
for jj=1:13
    v1 = dura_c_align(id_4g==1,jj)';
    v2 = dura_c_align(id_4g==2,jj)';
    v3 = dura_c_align(id_4g==3,jj)';
    v4 = dura_c_align(id_4g==4,jj)';
    
    val = [v1,v2,v3,v4]';
    group = [ones(1,length(v1))*1,ones(1,length(v2))*2,ones(1,length(v3))*3,ones(1,length(v4))*4]';
    ages = [br_list.Age(id_4g==1);br_list.Age(id_4g==2);br_list.Age(id_4g==3);br_list.Age(id_4g==4)];
    sex = [br_list.SEXE(id_4g==1);br_list.SEXE(id_4g==2);br_list.SEXE(id_4g==3);br_list.SEXE(id_4g==4)];
    
    tmptbl = table(val,categorical(group), ages,categorical(sex),'VariableNames',{'Val','Group','Age','Sex'});
    
    mdl = fitglm(tmptbl, 'Val ~ Group + Age + Sex');
    pval_c(1:3,jj) = mdl.Coefficients.pValue(2:4);
    
    tmptbl.Group = reordercats(tmptbl.Group,{'2','3','4','1'});
    mdl = fitglm(tmptbl, 'Val ~ Group + Age + Sex');
    pval_c(4:5,jj) = mdl.Coefficients.pValue(2:3);
    
    tmptbl.Group = reordercats(tmptbl.Group,{'3','4','1','2'});
    mdl = fitglm(tmptbl, 'Val ~ Group + Age + Sex');
    pval_c(6,jj) = mdl.Coefficients.pValue(2);
end

% Ajust

for ii=1:6
    pval_c_adj(ii,:) = mafdr(pval_c(ii,:),'BHFDR','true');
end

% Dynamic time warping of SWA over first main sleep episode

dtw_swa = nan(619,1000);

for ii=find(id_4g>0)'
   tab = spectral_tables{ii};
   if isempty(tab) | wake_up_index(ii)==0
       continue
   end
   
   hypno = tab.Hypno;
   time = tab.Time;
   cycles = tab.Cycles;
   wuid = wake_up_index(ii);
   id1 = find(cycles>0,1,'first');
   hypno_n1 = hypno(id1:wuid-1);
   cycles_n1 = cycles(id1:wuid-1);
   
   tmp = swa_array{ii}([3 5],:);
   tmp(:,tmp(2,:)>wuid/30/60) = [];
   
   swa_int = interp1((tmp(2,:)-tmp(2,1))/(tmp(2,end)-tmp(2,1)),tmp(1,:),linspace(0,1,1000),'linear');
    dtw_swa(ii,:) = swa_int;
end

f4c=figure
hold on
plot((1:1000)/10,nanmean(dtw_swa(id_4g==1,:)),'color',c_hi1,'LineWidth',2)
plot((1:1000)/10,nanmean(dtw_swa(id_4g==2,:)),'color',c_hi2,'LineWidth',2)
plot((1:1000)/10,nanmean(dtw_swa(id_4g==3,:)),'color',c_hi_nlst,'LineWidth',2)
plot((1:1000)/10,nanmean(dtw_swa(id_4g==4,:)),'color',c_tem,'LineWidth',2)
ylabel('Relative SWA')
title('Relative SWA as a function of N1')
xlabel('Proportion of the first main sleep episode (%)')


% Statistics : GLM per 10% window

pval_swa = [];
for ii=1:10
    tmp = mean(dtw_swa(:,(ii-1)*100+1:ii*100),2);
    h1 = lillietest(tmp(id_4g==1));
    h2 = lillietest(tmp(id_4g==2));
    h3 = lillietest(tmp(id_4g==3));
    h4 = lillietest(tmp(id_4g==4));
    
    tmptbl = table(tmp(id_4g>0), categorical(id_4g(id_4g>0)),br_list.Age(id_4g>0),categorical(br_list.SEXE(id_4g>0)),...
        'VariableNames',{'Val','Group','Age','Sex'});
    
    if h1+h2+h3+h4==0
        mdl = fitlm(tmptbl,'Val ~ Group + Age + Sex');
        pval_swa(1:3,ii) = mdl.Coefficients.pValue(2:4);
        
        tmptbl.Group = reordercats(tmptbl.Group,{'2','3','4','1'});
         mdl = fitlm(tmptbl,'Val ~ Group + Age + Sex');
        pval_swa(4:5,ii) = mdl.Coefficients.pValue(2:3);
        
        tmptbl.Group = reordercats(tmptbl.Group,{'3','4','1','2'});
         mdl = fitlm(tmptbl,'Val ~ Group + Age + Sex');
        pval_swa(6,ii) = mdl.Coefficients.pValue(2);
    else
        mdl = fitglm(tmptbl,'Val ~ Group + Age + Sex','Distribution','Gamma','link','log');
        pval_swa(1:3,ii) = mdl.Coefficients.pValue(2:4);
        
        tmptbl.Group = reordercats(tmptbl.Group,{'2','3','4','1'});
         mdl = fitglm(tmptbl,'Val ~ Group + Age + Sex','Distribution','Gamma','link','log');
        pval_swa(4:5,ii) = mdl.Coefficients.pValue(2:3);
        
        tmptbl.Group = reordercats(tmptbl.Group,{'3','4','1','2'});
         mdl = fitglm(tmptbl,'Val ~ Group + Age + Sex','Distribution','Gamma','link','log');
        pval_swa(6,ii) = mdl.Coefficients.pValue(2);
    end
end

% Adjust

pval_swa_adj = pval_swa;

pval_swa_adj(:) = mafdr(pval_swa_adj(:),'BHFDR','true')

 
%% FIGURE 5

% SWA rebound as a function of prior wake 

swa_rebounds_data = [];
for ii=find(id_4g>0)'
    tab = spectral_tables{ii};
    swa= tab.swa;
    cycles = tab.Cycles;
    id_wu = wake_up_index(ii);
    hypno = tab.Hypno;
    % Wake periods after main first sleep episode
    id_w = finddatagroups(cycles(id_wu:end),0) + id_wu;
    % Remove short sleeps (<45min)
    for jj = 2:2:length(id_w)-1
        if id_w(jj+1)-id_w(jj) < 45*30
            id_w(jj+1) = NaN;
            id_w(jj) = NaN;
        end
    end
    id_w(isnan(id_w)) = [];
    if id_w(end) > height(tab) - 45*30
        id_w(end) = NaN;
        id_w(end-1) = NaN;
    end
    id_w(isnan(id_w)) = [];
    
    % For each wake bout : SWA in following cycle
    for jj = 1:2:length(id_w)
        dura = (id_w(jj+1) - id_w(jj) +1)/30/60;
        next_c = cycles == cycles(id_w(jj+1)+1);
        swa_next_c = swa(next_c);
        hypno_next_c = hypno(next_c);
        swa_peak = prctile(swa_next_c(hypno_next_c>=4),95);      
        swa_rebounds_data = [swa_rebounds_data,[swa_peak;dura;ii]]; 
    end
end


xx = linspace(0,16,300);

f5=figure
hold on
scatter(swa_rebounds_data(2,ismember(swa_rebounds_data(3,:),find(id_4g==1))),swa_rebounds_data(1,ismember(swa_rebounds_data(3,:),find(id_4g==1))),[],c_hi1,'filled')
scatter(swa_rebounds_data(2,ismember(swa_rebounds_data(3,:),find(id_4g==2))),swa_rebounds_data(1,ismember(swa_rebounds_data(3,:),find(id_4g==2))),[],c_hi2,'filled')
scatter(swa_rebounds_data(2,ismember(swa_rebounds_data(3,:),find(id_4g==3))),swa_rebounds_data(1,ismember(swa_rebounds_data(3,:),find(id_4g==3))),[],c_hi_nlst,'filled')
scatter(swa_rebounds_data(2,ismember(swa_rebounds_data(3,:),find(id_4g==4))),swa_rebounds_data(1,ismember(swa_rebounds_data(3,:),find(id_4g==4))),[],c_tem,'filled')
xlabel('Prior wake duration (hours)')
ylabel('Relative SWA')
title('95th percentile of SWA in the sleep cycle following wake')
% Fit linear trends
% 1
x1=swa_rebounds_data(2,ismember(swa_rebounds_data(3,:),find(id_4g==1)))';
y1=swa_rebounds_data(1,ismember(swa_rebounds_data(3,:),find(id_4g==1)))';
tmpage = br_list.Age(swa_rebounds_data(3,ismember(swa_rebounds_data(3,:),find(id_4g==1))));
tmpsex = categorical(br_list.SEXE(swa_rebounds_data(3,ismember(swa_rebounds_data(3,:),find(id_4g==1)))));
idnan = isnan(x1) | isnan(y1);
x1(idnan) = [];
y1(idnan) = [];
tmpage(idnan) = [];
tmpsex(idnan) = [];
tbl = table(x1,y1,tmpage,tmpsex,...
    'VariableNames',{'X','Y','Age','Sex'});
mdl1=fitlm(tbl, 'Y ~X + Age + Sex');
plot(xx,mdl1.Coefficients.Estimate(1) + xx * mdl1.Coefficients.Estimate(2),'color',c_hi1,'LineWidth',2)
% 2
x2=swa_rebounds_data(2,ismember(swa_rebounds_data(3,:),find(id_4g==2)))';
y2=swa_rebounds_data(1,ismember(swa_rebounds_data(3,:),find(id_4g==2)))';
tmpage = br_list.Age(swa_rebounds_data(3,ismember(swa_rebounds_data(3,:),find(id_4g==2))));
tmpsex = categorical(br_list.SEXE(swa_rebounds_data(3,ismember(swa_rebounds_data(3,:),find(id_4g==2)))));
idnan = isnan(x2) | isnan(y2);
x2(idnan) = [];
y2(idnan) = [];
tmpage(idnan) = [];
tmpsex(idnan) = [];
tbl = table(x2,y2,tmpage,tmpsex,...
    'VariableNames',{'X','Y','Age','Sex'});
mdl2=fitlm(tbl, 'Y ~X + Age + Sex');
plot(xx,mdl2.Coefficients.Estimate(1) + xx * mdl2.Coefficients.Estimate(2),'color',c_hi2,'LineWidth',2)
% 3
x3=swa_rebounds_data(2,ismember(swa_rebounds_data(3,:),find(id_4g==3)))';
y3=swa_rebounds_data(1,ismember(swa_rebounds_data(3,:),find(id_4g==3)))';
tmpage = br_list.Age(swa_rebounds_data(3,ismember(swa_rebounds_data(3,:),find(id_4g==3))));
tmpsex = categorical(br_list.SEXE(swa_rebounds_data(3,ismember(swa_rebounds_data(3,:),find(id_4g==3)))));
idnan = isnan(x3) | isnan(y3);
x3(idnan) = [];
y3(idnan) = [];
tmpage(idnan) = [];
tmpsex(idnan) = [];
tbl = table(x3,y3,tmpage,tmpsex,...
    'VariableNames',{'X','Y','Age','Sex'});
mdl3=fitlm(tbl, 'Y ~X + Age + Sex')
plot(xx,mdl3.Coefficients.Estimate(1) + xx * mdl3.Coefficients.Estimate(2),'color',c_hi_nlst,'LineWidth',2)
% 4
x4=swa_rebounds_data(2,ismember(swa_rebounds_data(3,:),find(id_4g==4)))';
y4=swa_rebounds_data(1,ismember(swa_rebounds_data(3,:),find(id_4g==4)))';
tmpage = br_list.Age(swa_rebounds_data(3,ismember(swa_rebounds_data(3,:),find(id_4g==4))));
tmpsex = categorical(br_list.SEXE(swa_rebounds_data(3,ismember(swa_rebounds_data(3,:),find(id_4g==4)))));
idnan = isnan(x4) | isnan(y4);
x4(idnan) = [];
y4(idnan) = [];
tmpage(idnan) = [];
tmpsex(idnan) = [];
tbl = table(x4,y4,tmpage,tmpsex,...
    'VariableNames',{'X','Y','Age','Sex'});
mdl4=fitlm(tbl, 'Y ~X + Age + Sex')
plot(xx,mdl4.Coefficients.Estimate(1) + xx * mdl4.Coefficients.Estimate(2),'color',c_tem,'LineWidth',2)


tmpgroup = ismember(swa_rebounds_data(3,:),find(id_4g==1))*1+...
    ismember(swa_rebounds_data(3,:),find(id_4g==2))*2+...
    ismember(swa_rebounds_data(3,:),find(id_4g==3))*3+...
    ismember(swa_rebounds_data(3,:),find(id_4g==4))*4;

% STATS : GLME 
tmptbl = array2table([swa_rebounds_data',tmpgroup'],'VariableNames',{'SWA','duration','Subject','Group'});
tmptbl.Age = br_list.Age(tmptbl.Subject);
tmptbl.Sex = categorical(br_list.SEXE(tmptbl.Subject));
tmptbl.Subject = categorical(tmptbl.Subject);
tmptbl.Group = categorical(tmptbl.Group);

pval = [];
mdl = fitglme(tmptbl, 'SWA ~ duration*Group + Age + Sex +(1|Subject)')
pval(1:3) = mdl.Coefficients.pValue(8:10);
tmptbl.Group = reordercats(tmptbl.Group,{'2','3','4','1'});
mdl = fitglme(tmptbl, 'SWA ~ duration*Group + Age + Sex +(1|Subject)')
pval(4:5) = mdl.Coefficients.pValue(8:9);
tmptbl.Group = reordercats(tmptbl.Group,{'3','4','1','2'});
mdl = fitglme(tmptbl, 'SWA ~ duration*Group + Age + Sex +(1|Subject)')
pval(6) = mdl.Coefficients.pValue(8);
pval_adj = mafdr(pval,'BHFDR','true')
tmptbl.Group = reordercats(tmptbl.Group,{'4','1','2','3'});
mdl = fitglme(tmptbl, 'SWA ~ duration*Group + Age + Sex +(1|Subject)')
% Difference in slope 1,2>3,4

%% FIT MATHEMATICAL MODEL

% Reminder, description of the model is:
% S(0) = So
% S(t+dt) = S(t) + [UA - (UA-S(t)) * exp(-dt / tr_stage)] - [(S(t)-LA) * exp(-dt/td_stage) + LA]
% where: 
% S(t) is the sleep pressure across time
% UA = upper asypmtote
% LA = lower asymptote
% tr_stage is the stage-dependent rise constant
% td_stage is the stage-dependent decay constant

sol_modl1_v3 = {};

for ii=find(id_4g>0)'
    ii
    tab = spectral_tables{ii};
    hypno = tab.Hypno;
    idfirst = find(hypno>1,1,'first');
    tab = tab(idfirst:end,:);
    swa = tab.swa;
    cycles = tab.Cycles;
    hypno = tab.Hypno;
    time = tab.Time;
    time_h = hours(time);
    % Get SWA at each midcycle: 3 options : median SWA over cycle, 95th
    % percentile, or AUC
    %swa_array: 1 = swa median, 2 = swa auc, 3 = sa max, 4 = time midcycle, 5 = time start cycle
    % Here we select 95th percentile and start of cycle
    obs_P = swa_array{ii}(3,:);
    obs_t = swa_array{ii}(5,:) - idfirst/30/60/60;
    obs_t_mid = swa_array{ii}(4,:) - idfirst/30/60/60;
    obs_t(isnan(obs_P)) = [];
    obs_P(isnan(obs_P)) = [];
    
    swa_sws = swa(hypno>3);
   

    % First : remove 1st SWA points if lower than second
    
    if obs_P(1)<obs_P(2)
       obs_P(1) = [];
       obs_t(1) = [];
    end
    
    % Then: remove sleep bouts that are shorter than 90min and dont reach
    % N3
    
    sb = finddatagroups(cycles>0,1);
    for jj=1:2:length(sb)
        idstart = sb(jj);
        idend = sb(jj+1);
        hypno_sb = hypno(idstart:idend);
        time_sb = time_h(idstart:idend);
        if idend - idstart +1 < 90*30 & sum(hypno_sb==5)==0
            id_rm = obs_t>=time_sb(1) & obs_t<=time_sb(end);
            obs_t(id_rm) = [];
            obs_P(id_rm) = [];
        end
    end
    
    
% MODEL  
% function that solves the model was coded in fit_S_dual_stages.m file
% which is available on the github as well
swa_sp = swa(hypno==2);
% Set initial values (1st column) and bounds (2nd lower, 3rd higher)
hyperparameters = [1, prctile(swa_sws,1),max(swa_sws);... % SO
    4, 0.1, 50;...   % Td_nrem
    50, 0.1, 50;...   % Tr_nrem   
    4, 0.1, 50;...   % Td_rem
    20, 0.1, 50;...   % Tr_rem
    20, 0.1, 50;...   % Td_w
    5, 0.1, 50;...   % Tr_w
    0.95*min(swa_sp), 0.95*min(swa_sp), 0.95*min(swa_sp);... %LA
    prctile(swa_sws,99.9),prctile(swa_sws,99.9),1.2*prctile(swa_sws,99.9)]; %UA
% Run optimizer fmincon
options = optimoptions('fmincon');
% Set linear constraints
A = [0,1,0,-1,0,0,0,0,0;0,0,-1,0,1,0,0,0,0;0,0,0,0,0,-1,1,0,0];
b = [0;0;0];
[xval,fval] = fmincon(@(params) fit_S_dual_stages(params,hypno, time_h, obs_t, obs_P),...
    hyperparameters(:,1),...
    A,b,[],[],hyperparameters(:,2),hyperparameters(:,3),[],...
    options);


% Obtain optimized sleep pressure timecourse
[~,val_f_mdl2] = fit_S_dual_stages(xval, hypno, time_h, obs_t, obs_P);

% Store results in structure
parameters = struct('SO',xval(1),'Td_nrem',xval(2),'Tr_nrem',xval(3),'Td_rem',xval(4),'Tr_rem',xval(5),'Td_w',xval(6),'Tr_w',xval(7),'LA',xval(8),'UA',xval(9));
sol = struct('parameters',parameters,...
    'error',fval,...
    'hypno',hypno,...
    'func',val_f_mdl2);

sol_modl1_v3{ii} = sol;
   
end


all_errs3 = nan(1,619);
all_errs3(id_4g>0) = cellfun(@(x) x.error,sol_modl1_v3(id_4g>0));
% To visualize distribution and set the threshold on fit error; we selected
% 1.2

% Stats on parameters

ptot = [];
% SO
val = [cellfun(@(x) x.parameters.SO,sol_modl1_v3(id_4g==1 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.SO,sol_modl1_v3(id_4g==2 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.SO,sol_modl1_v3(id_4g==3 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.SO,sol_modl1_v3(id_4g==4 & all_errs3'<=1.2))]';
groups = [ones(1,sum(id_4g==1 & all_errs3'<=1.2))*1,...
    ones(1,sum(id_4g==2 & all_errs3'<=1.2))*2,...
    ones(1,sum(id_4g==3 & all_errs3'<=1.2))*3,...
    ones(1,sum(id_4g==4 & all_errs3'<=1.2))*4]';

age = [br_list.Age(id_4g==1 & all_errs3'<=1.2);...
    br_list.Age(id_4g==2 & all_errs3'<=1.2);...
    br_list.Age(id_4g==3 & all_errs3'<=1.2);...
    br_list.Age(id_4g==4 & all_errs3'<=1.2)];
sex = categorical([br_list.SEXE(id_4g==1 & all_errs3'<=1.2);...
    br_list.SEXE(id_4g==2 & all_errs3'<=1.2);...
    br_list.SEXE(id_4g==3 & all_errs3'<=1.2);...
    br_list.SEXE(id_4g==4 & all_errs3'<=1.2)]);

tabtmp = table(val,groups,age,sex,'VariableNames',{'Value','Group','Age','Sex'});
tabtmp.Group = categorical(tabtmp.Group);
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Gamma','Link', 'log');
mdl_reduced = fitglme(tabtmp, 'Value ~ Age + Sex', 'Distribution', 'Gamma', 'Link', 'log');
compare(mdl,mdl_reduced)
ptot(1:3,1) = mdl.Coefficients.pValue(2:4);
tabtmp.Group = reordercats(tabtmp.Group,{'2','3','4','1'});
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Gamma','Link', 'log');
ptot(4:5,1) = mdl.Coefficients.pValue(2:3);
tabtmp.Group = reordercats(tabtmp.Group,{'3','4','1','2'});
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Gamma','Link', 'log');
ptot(6,1) = mdl.Coefficients.pValue(2);

ptot(:,1) = mafdr(ptot(:,1),'BHFDR','true')

% UA
val = [cellfun(@(x) x.parameters.UA,sol_modl1_v3(id_4g==1 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.UA,sol_modl1_v3(id_4g==2 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.UA,sol_modl1_v3(id_4g==3 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.UA,sol_modl1_v3(id_4g==4 & all_errs3'<=1.2))]';

tabtmp = table(val,groups,age,sex,'VariableNames',{'Value','Group','Age','Sex'});
tabtmp.Group = categorical(tabtmp.Group);
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Gamma','Link', 'log');
mdl_reduced = fitglme(tabtmp, 'Value ~ Age + Sex', 'Distribution', 'Gamma', 'Link', 'log');
compare(mdl,mdl_reduced)
ptot(1:3,2) = mdl.Coefficients.pValue(2:4);
tabtmp.Group = reordercats(tabtmp.Group,{'2','3','4','1'});
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Gamma','Link', 'log');
ptot(4:5,2) = mdl.Coefficients.pValue(2:3);
tabtmp.Group = reordercats(tabtmp.Group,{'3','4','1','2'});
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Gamma','Link', 'log');
ptot(6,2) = mdl.Coefficients.pValue(2);
ptot(:,2) = mafdr(ptot(:,2),'BHFDR','true')

% LA
val = [cellfun(@(x) x.parameters.LA,sol_modl1_v3(id_4g==1 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.LA,sol_modl1_v3(id_4g==2 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.LA,sol_modl1_v3(id_4g==3 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.LA,sol_modl1_v3(id_4g==4 & all_errs3'<=1.2))]';

tabtmp = table(val,groups,age,sex,'VariableNames',{'Value','Group','Age','Sex'});
tabtmp.Group = categorical(tabtmp.Group);
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Normal');
mdl_reduced = fitglme(tabtmp, 'Value ~ Age + Sex', 'Distribution', 'Normal');
compare(mdl,mdl_reduced)
ptot(1:3,3) = mdl.Coefficients.pValue(2:4);
tabtmp.Group = reordercats(tabtmp.Group,{'2','3','4','1'});
mdl = fitglme(tabtmp,'Value ~ Group+ Age + Sex','Distribution','Normal');
ptot(4:5,3) = mdl.Coefficients.pValue(2:3);
tabtmp.Group = reordercats(tabtmp.Group,{'3','4','1','2'});
mdl = fitglme(tabtmp,'Value ~ Group+ Age + Sex','Distribution','Normal');
ptot(6,3) = mdl.Coefficients.pValue(2);
ptot(:,3) = mafdr(ptot(:,3),'BHFDR','true')

% Td NREM
val = [cellfun(@(x) x.parameters.Td_nrem,sol_modl1_v3(id_4g==1 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Td_nrem,sol_modl1_v3(id_4g==2 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Td_nrem,sol_modl1_v3(id_4g==3 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Td_nrem,sol_modl1_v3(id_4g==4 & all_errs3'<=1.2))]';

[p,~,stats] = kruskalwallis(val,groups);
c=multcompare(stats)
ptot(:,4) = c(:,6);

% Tr NREM
val = [cellfun(@(x) x.parameters.Tr_nrem,sol_modl1_v3(id_4g==1 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Tr_nrem,sol_modl1_v3(id_4g==2 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Tr_nrem,sol_modl1_v3(id_4g==3 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Tr_nrem,sol_modl1_v3(id_4g==4 & all_errs3'<=1.2))]';
[p,~,stats] = kruskalwallis(val,groups);
c=multcompare(stats)
ptot(:,5) = c(:,6);

% Td REM

val = [cellfun(@(x) x.parameters.Td_rem,sol_modl1_v3(id_4g==1 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Td_rem,sol_modl1_v3(id_4g==2 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Td_rem,sol_modl1_v3(id_4g==3 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Td_rem,sol_modl1_v3(id_4g==4 & all_errs3'<=1.2))]';
[p,~,stats] = kruskalwallis(val,groups);
c=multcompare(stats)
ptot(:,6) = c(:,6);

% Tr REM

val = [cellfun(@(x) x.parameters.Tr_rem,sol_modl1_v3(id_4g==1 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Tr_rem,sol_modl1_v3(id_4g==2 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Tr_rem,sol_modl1_v3(id_4g==3 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Tr_rem,sol_modl1_v3(id_4g==4 & all_errs3'<=1.2))]';
[p,~,stats] = kruskalwallis(val,groups);
c=multcompare(stats)
ptot(:,7) = c(:,6);

% Td wake

val = [cellfun(@(x) x.parameters.Td_w,sol_modl1_v3(id_4g==1 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Td_w,sol_modl1_v3(id_4g==2 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Td_w,sol_modl1_v3(id_4g==3 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Td_w,sol_modl1_v3(id_4g==4 & all_errs3'<=1.2))]';

[mean(val(groups==1)) std(val(groups==1));...
    mean(val(groups==2)) std(val(groups==2));...
    mean(val(groups==3)) std(val(groups==3));...
    mean(val(groups==4)) std(val(groups==4))]

tabtmp = table(val>=49.5,groups,age,sex,'VariableNames',{'Value','Group','Age','Sex'});
tabtmp.Group = categorical(tabtmp.Group);
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Binomial');
ptot(1:3,8) = mdl.Coefficients.pValue(2:4);
tabtmp.Group = reordercats(tabtmp.Group,{'2','3','4','1'});
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Binomial');
ptot(4:5,8) = mdl.Coefficients.pValue(2:3);
tabtmp.Group = reordercats(tabtmp.Group,{'3','4','1','2'});
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Binomial');
ptot(6,8) = mdl.Coefficients.pValue(2);
ptot(:,8) = mafdr(ptot(:,8),'BHFDR','true')

[p,~,stats] = kruskalwallis(val,groups);
c=multcompare(stats)
ptot(:,8) = c(:,6);
% Tr wake

val = [cellfun(@(x) x.parameters.Tr_w,sol_modl1_v3(id_4g==1 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Tr_w,sol_modl1_v3(id_4g==2 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Tr_w,sol_modl1_v3(id_4g==3 & all_errs3'<=1.2)),...
    cellfun(@(x) x.parameters.Tr_w,sol_modl1_v3(id_4g==4 & all_errs3'<=1.2))]';

[mean(val(groups==1)) std(val(groups==1));...
    mean(val(groups==2)) std(val(groups==2));...
    mean(val(groups==3)) std(val(groups==3));...
    mean(val(groups==4)) std(val(groups==4))]

tabtmp = table(val>=10,groups,age,sex,'VariableNames',{'Value','Group','Age','Sex'});
tabtmp.Group = categorical(tabtmp.Group);
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Binomial');
ptot(1:3,9) = mdl.Coefficients.pValue(2:4);
tabtmp.Group = reordercats(tabtmp.Group,{'2','3','4','1'});
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Binomial');
ptot(4:5,9) = mdl.Coefficients.pValue(2:3);
tabtmp.Group = reordercats(tabtmp.Group,{'3','4','1','2'});
mdl = fitglme(tabtmp,'Value ~ Group + Age + Sex','Distribution','Binomial');
ptot(6,9) = mdl.Coefficients.pValue(2);
ptot(:,9) = mafdr(ptot(:,9),'BHFDR','true')

[p,~,stats] = kruskalwallis(val,groups);
c=multcompare(stats)
ptot(:,9) = c(:,6);


%% FIGURE 6 

% Here we investigate the angle of change of sleep pressure upon 2min in
% the different sleep stages, at different starting pressure levels

dur = 2;
all_rem = ones(1,dur*30,11)*2;
all_nrem = ones(1,dur*30,11)*3;
all_w = ones(1,dur*30,11)*1;
time_h = (1:dur*30)/30/60;
mat_all_nrem = nan(619,dur*30,11);
mat_all_rem = nan(619,dur*30,11);
mat_all_w = nan(619,dur*30,11);

slope_nrem = nan(619,11);
slope_rem = nan(619,11);
slope_w = nan(619,11);

k=1
for jj=[1 10 20 30 40 50 60 70 80 90 99]/100 % Pressure levels evaluated
    for ii=find(id_4g>0 & all_errs3'<=1.2)'
        sol = sol_modl1_v3{ii};
        params_sleep = [sol.parameters.UA*jj,sol.parameters.Td_nrem,sol.parameters.Tr_nrem,sol.parameters.Td_rem,sol.parameters.Tr_rem,sol.parameters.Td_w,sol.parameters.Tr_w,sol.parameters.LA,sol.parameters.UA];
        params_w = [sol.parameters.UA*jj,sol.parameters.Td_nrem,sol.parameters.Tr_nrem,sol.parameters.Td_rem,sol.parameters.Tr_rem,sol.parameters.Td_w,sol.parameters.Tr_w,sol.parameters.LA,sol.parameters.UA];
        [~,sim_S_nrem] = fit_S_dual_stages(params_sleep, all_nrem, time_h);
        [~,sim_S_rem] = fit_S_dual_stages(params_sleep, all_rem, time_h);
        [~,sim_S_w] = fit_S_dual_stages(params_w, all_w, time_h);
        
        mat_all_nrem(ii,:,k) = sim_S_nrem;
        mat_all_rem(ii,:,k) = sim_S_rem;
        mat_all_w(ii,:,k) = sim_S_w;
        
        mdl = fitlm(time_h,sim_S_nrem);
        slope_nrem(ii,k) = mdl.Coefficients.Estimate(2);
        mdl = fitlm(time_h,sim_S_rem);
        slope_rem(ii,k) = mdl.Coefficients.Estimate(2);
        mdl = fitlm(time_h,sim_S_w);
        slope_w(ii,k) = mdl.Coefficients.Estimate(2);
        
    end
    k=k+1;
end

% Version compass plot
% From each slope, compute angle, then get average angle for each group

theta_nrem = atan(slope_nrem);
theta_rem = atan(slope_rem);
theta_w = atan(slope_w);

f6=figure
t=tiledlayout(3,6,'TileSpacing','none')
pressures = [1 10 20 30 40 50 60 70 80 90 99];
tmp = [1 3:2:11]; % Only plot some pressure levels for legibility
for i=1:6
    jj=tmp(i);
    mm = tmp(7-i);
    t1=nexttile(i);
    [x1,y1]=pol2cart(nanmean(theta_nrem(id_4g==1 & all_errs3'<=1.2,mm)),1);
    [x2,y2]=pol2cart(nanmean(theta_nrem(id_4g==2 & all_errs3'<=1.2,mm)),1);
    [x3,y3]=pol2cart(nanmean(theta_nrem(id_4g==3 & all_errs3'<=1.2,mm)),1);
    [x4,y4]=pol2cart(nanmean(theta_nrem(id_4g==4 & all_errs3'<=1.2,mm)),1);    
    c=compass([x1,x2,x3,x4],[y1,y2,y3,y4]);
    c(1).Color = c_hi1;
    c(2).Color = c_hi2;
    c(3).Color = c_hi_nlst;
    c(4).Color = c_tem;
    for kk=1:4;c(kk).LineWidth=1.0;end;
    xlim([0 1])
    h = findall(gca,'Type','Text');for ii=1:length(h);h(ii).String = '';end
    title(t1,strcat(string(pressures(mm)),'%'))
    if i==1;ylabel('Angle of change - NREM');end
    
    nexttile(i+6);
    [x1,y1]=pol2cart(nanmean(theta_rem(id_4g==1 & all_errs3'<=1.2,mm)),1);
    [x2,y2]=pol2cart(nanmean(theta_rem(id_4g==2 & all_errs3'<=1.2,mm)),1);
    [x3,y3]=pol2cart(nanmean(theta_rem(id_4g==3 & all_errs3'<=1.2,mm)),1);
    [x4,y4]=pol2cart(nanmean(theta_rem(id_4g==4 & all_errs3'<=1.2,mm)),1);    
    c=compass([x1,x2,x3,x4],[y1,y2,y3,y4]);
    c(1).Color = c_hi1;
    c(2).Color = c_hi2;
    c(3).Color = c_hi_nlst;
    c(4).Color = c_tem;
    for kk=1:4;c(kk).LineWidth=1.0;end;
    xlim([0 1])
     h = findall(gca,'Type','Text');for ii=1:length(h);h(ii).String = '';end
     if i==1;ylabel('Angle of change - REM');end
     
    t3=nexttile(i+12);
    [x1,y1]=pol2cart(nanmean(theta_w(id_4g==1 & all_errs3'<=1.2,jj)),1);
    [x2,y2]=pol2cart(nanmean(theta_w(id_4g==2 & all_errs3'<=1.2,jj)),1);
    [x3,y3]=pol2cart(nanmean(theta_w(id_4g==3 & all_errs3'<=1.2,jj)),1);
    [x4,y4]=pol2cart(nanmean(theta_w(id_4g==4 & all_errs3'<=1.2,jj)),1);    
    c=compass([x1,x2,x3,x4],[y1,y2,y3,y4]);
    c(1).Color = c_hi1;
    c(2).Color = c_hi2;
    c(3).Color = c_hi_nlst;
    c(4).Color = c_tem;
    for kk=1:4;c(kk).LineWidth=1.0;end;
    xlim([0 1])
     h = findall(gca,'Type','Text');for ii=1:length(h);h(ii).String = '';end
     title(t3,strcat(string(pressures(jj)),'%'))
      if i==1;ylabel('Angle of change - wake');end
end
title(t,{'Direction of the dynamics in the different vigilance stages as a function of sleep pressure',''})

% STATS - GLME

pval_nrem2 = [];
pressures_tbl = repelem(pressures,215,1);
tmp = theta_nrem(id_4g>0 & all_errs3'<=1.2,:);
groups_tbl = repelem(id_4g(id_4g>0 & all_errs3'<=1.2),1,11);
ages = repelem(br_list.Age(id_4g>0 & all_errs3'<=1.2),1,11);
sexes = repelem(br_list.SEXE(id_4g>0 & all_errs3'<=1.2),1,11);
tbltmp = table(tmp(:),groups_tbl(:),pressures_tbl(:),ages(:),categorical(sexes(:)),...
    'VariableNames',{'Theta','Group','Pressure','Age','Sex'});
tbltmp.Group = categorical(tbltmp.Group);
mdl_nrem = fitglme(tbltmp,'Theta~Group*Pressure + Age + Sex','Distribution','Normal');
pval_nrem2(1:3) = mdl_nrem.Coefficients.pValue(8:10);
tbltmp.Group = reordercats(tbltmp.Group,{'2','3','4','1'});
mdl_nrem = fitglme(tbltmp,'Theta~Group*Pressure + Age + Sex','Distribution','Normal');
pval_nrem2(4:5) = mdl_nrem.Coefficients.pValue(8:9);
tbltmp.Group = reordercats(tbltmp.Group,{'3','4','1','2'});
mdl_nrem = fitglme(tbltmp,'Theta~Group*Pressure + Age + Sex','Distribution','Normal');
pval_nrem2(6) = mdl_nrem.Coefficients.pValue(8);

pval_rem2 = [];
tmp = theta_rem(id_4g>0 & all_errs3'<=1.2,:);
tbltmp = table(tmp(:),groups_tbl(:),pressures_tbl(:),ages(:),categorical(sexes(:)),...
    'VariableNames',{'Theta','Group','Pressure','Age','Sex'});
tbltmp.Group = categorical(tbltmp.Group);
mdl_rem = fitglme(tbltmp,'Theta~Group*Pressure + Age + Sex','Distribution','Normal');
pval_rem2(1:3) = mdl_rem.Coefficients.pValue(8:10);
tbltmp.Group = reordercats(tbltmp.Group,{'2','3','4','1'});
mdl_rem = fitglme(tbltmp,'Theta~Group*Pressure + Age + Sex','Distribution','Normal');
pval_rem2(4:5) = mdl_rem.Coefficients.pValue(8:9);
tbltmp.Group = reordercats(tbltmp.Group,{'3','4','1','2'});
mdl_rem = fitglme(tbltmp,'Theta~Group*Pressure + Age + Sex','Distribution','Normal');
pval_rem2(6) = mdl_rem.Coefficients.pValue(8);

pval_w2 = [];
tmp = theta_w(id_4g>0 & all_errs3'<=1.2,:);
tbltmp = table(tmp(:),groups_tbl(:),pressures_tbl(:),ages(:),categorical(sexes(:)),...
    'VariableNames',{'Theta','Group','Pressure','Age','Sex'});
tbltmp.Group = categorical(tbltmp.Group);
mdl_w = fitglme(tbltmp,'Theta~Group*Pressure + Age + Sex','Distribution','Normal');
pval_w2(1:3) = mdl_w.Coefficients.pValue(8:10);
tbltmp.Group = reordercats(tbltmp.Group,{'2','3','4','1'});
mdl_w = fitglme(tbltmp,'Theta~Group*Pressure + Age + Sex','Distribution','Normal');
pval_w2(4:5) = mdl_w.Coefficients.pValue(8:9);
tbltmp.Group = reordercats(tbltmp.Group,{'3','4','1','2'});
mdl_w = fitglme(tbltmp,'Theta~Group*Pressure + Age + Sex','Distribution','Normal');
pval_w2(6) = mdl_w.Coefficients.pValue(8);


%% FIGURE 7

% Next we want to simulate sleep pressure to standardize sleep-wake patterns across gorups
% But we want to keep the NREM/REM proportions of each individual: extract that information 

prop_rem = nan(1,619);
prop_w_n1 = nan(1,619);

for ii=find(id_4g>0)'
   tab = spectral_tables{ii};
   if isempty(tab)
       continue
   end
    wuid = wake_up_index(ii);
    if wuid==0
        continue
    end
    hypno = tab.Hypno(1:wu_id);
    prop_rem = sum(hypno==2) / sum(hypno>1) * 100;
    prop_rem(ii) = prop_rem;
    prop_w_n1(ii) = sum(hypno==1) / sum(hypno>1) * 100;
end

% Evaluate S in extended sleep

S_extended_sleep3 = [];
S_extended_sleep_boxplot3 = [];
time_h = (1:length(sw_extended))/30/60;
for ii=find(id_4g>0 & (all_errs3<=1.2)')'
    ii
    prop = prop_rem(ii);
    sol = sol_modl1_v3{ii};
    params = [sol.parameters.SO,sol.parameters.Td_nrem,sol.parameters.Tr_nrem,sol.parameters.Td_rem,sol.parameters.Tr_rem,sol.parameters.Td_w,sol.parameters.Tr_w,sol.parameters.LA,sol.parameters.UA];
    sw_extended = generateHypno(12*60*30,prop); % To generate the simulated NREM/REM timecourse while keeping individual REM proportion
    [~,val_S] = fit_S_dual_stages(params, sw_extended, time_h);
    % Extract value at each cycle start
    S_extended_sleep3(ii,:) = val_S;
    S_extended_sleep_boxplot3 = [S_extended_sleep_boxplot3, [[val_S(1:90*30:end);val_S(end)]';ones(1,9)*id_4g(ii);0:1.5:12]];
end

S_extended_sleep_boxplot3(3,S_extended_sleep_boxplot3(2,:)==2) = S_extended_sleep_boxplot3(3,S_extended_sleep_boxplot3(2,:)==2)+0.2;
S_extended_sleep_boxplot3(3,S_extended_sleep_boxplot3(2,:)==3) = S_extended_sleep_boxplot3(3,S_extended_sleep_boxplot3(2,:)==3)+0.4;
S_extended_sleep_boxplot3(3,S_extended_sleep_boxplot3(2,:)==4) = S_extended_sleep_boxplot3(3,S_extended_sleep_boxplot3(2,:)==4)+0.6;

pos = repelem(0:8,1,4);
pos(2:4:end) = pos(2:4:end)+0.2;
pos(3:4:end) = pos(3:4:end)+0.4;
pos(4:4:end) = pos(4:4:end)+0.6;

f7a=figure
hold on
boxplot(S_extended_sleep_boxplot3(1,:),S_extended_sleep_boxplot3(3,:),'Symbol','','Positions',pos)
plot(time_h,mean(S_extended_sleep3(id_4g==1 & all_errs3'<=1.2,:),1),'color',c_hi1,'LineWidth',2)
plot(time_h,mean(S_extended_sleep3(id_4g==2 & all_errs3'<=1.2,:),1),'color',c_hi2,'LineWidth',2)
plot(time_h,mean(S_extended_sleep3(id_4g==3 & all_errs3'<=1.2,:),1),'color',c_hi_nlst,'LineWidth',2)
plot(time_h,mean(S_extended_sleep3(id_4g==4 & all_errs3'<=1.2,:),1),'color',c_tem,'LineWidth',2)
h = findobj(gca,'Tag','Box')
for ii=1:4:36; p4=patch(get(h(ii),'XData'),get(h(ii),'YData'),c_tem,'FaceAlpha',.7);end
for ii=2:4:36; p3=patch(get(h(ii),'XData'),get(h(ii),'YData'),c_hi_nlst,'FaceAlpha',.7);end
for ii=3:4:36; p2=patch(get(h(ii),'XData'),get(h(ii),'YData'),c_hi2,'FaceAlpha',.7);end
for ii=4:4:36; p1=patch(get(h(ii),'XData'),get(h(ii),'YData'),c_hi1,'FaceAlpha',.7);end
xticks((0:9)+0.3)
xticklabels(string(0:1.5:12))
xlabel('Sleep duration (hours)')
ylabel('Theoretical sleep pressure')
legend([p1 p2 p3 p4],{'IH-LST-S','IH-LST-nS','IH-nLST','HC'})
title('Simulated sleep pressure during extended sleep')
grid on

% STATS : GLME per hour

pvalues_extended_glme3 = [];
k=1;
for ii=0:1.5:12
    val = S_extended_sleep_boxplot3(1,S_extended_sleep_boxplot3(3,:)-mod(S_extended_sleep_boxplot3(3,:),1.5)==ii);
    group = S_extended_sleep_boxplot3(2,S_extended_sleep_boxplot3(3,:)-mod(S_extended_sleep_boxplot3(3,:),1.5)==ii);
    tts = br_list.BEDREST_TTS32_h(id_4g>0 & all_errs3'<=1.2);
    age = br_list.Age(id_4g>0 & all_errs3'<=1.2);
    sex = br_list.SEXE(id_4g>0 & all_errs3'<=1.2);
    wuid = wake_up_index(id_4g>0 & all_errs3'<=1.2)/30/60;
    tmptbl = table(val',group',wuid,age,categorical(sex),'VariableNames',{'S','Group','TTS_N1','Age','Sex'});
    tmptbl.Group = categorical(tmptbl.Group);
    mdl = fitglme(tmptbl, 'S ~ Group + TTS_N1 + Age + Sex');
    pvalues_extended_glme3(1:3,k) = mdl.Coefficients.pValue(2:4);
    tmptbl.Group = reordercats(tmptbl.Group,{'2','3','4','1'});
    mdl = fitglme(tmptbl, 'S ~ Group + TTS_N1 + Age + Sex');
    pvalues_extended_glme3(4:5,k) = mdl.Coefficients.pValue(2:3);
    tmptbl.Group = reordercats(tmptbl.Group,{'3','4','1','2'});
    mdl = fitglme(tmptbl, 'S ~ Group + TTS_N1 + Age + Sex');
    pvalues_extended_glme3(6,k) = mdl.Coefficients.pValue(2);
    
    % Adjust
    pvalues_extended_glme3(:,k) = mafdr(pvalues_extended_glme3(:,k),'BHFDR',true);
    k=k+1;
end


% Simulate S during wake, after 8 hours of sleep

S_post_wake3 = [];
S_post_wake_boxplot3 = [];
for ii=find(id_4g>0 & (all_errs3<=1.2)')'
    ii
    prop = prop_rem_n1(ii);
    sol = sol_modl1_v3{ii};
    params = [sol.parameters.SO,sol.parameters.Td_nrem,sol.parameters.Tr_nrem,sol.parameters.Td_rem,sol.parameters.Tr_rem,sol.parameters.Td_w,sol.parameters.Tr_w,sol.parameters.LA,sol.parameters.UA];
    sw_extended = generateHypno(12*60*30,prop);
    sw_post_wake = [sw_extended(1:8*60*30) ones(1,16*60*30)];
    time_h = (1:length(sw_post_wake))/30/60;
    [~,val_S] = fit_S_dual_stages(params, sw_post_wake, time_h);
    S_post_wake3(ii,:) = val_S(9*30*60:end);
    S_post_wake_boxplot3 = [S_post_wake_boxplot3, [val_S(9*30*60:30*60:end)';ones(1,16)*id_4g(ii);1:16]];
end
S_post_wake_boxplot3(3,S_post_wake_boxplot3(2,:)==2) = S_post_wake_boxplot3(3,S_post_wake_boxplot3(2,:)==2)+0.2;
S_post_wake_boxplot3(3,S_post_wake_boxplot3(2,:)==3) = S_post_wake_boxplot3(3,S_post_wake_boxplot3(2,:)==3)+0.4;
S_post_wake_boxplot3(3,S_post_wake_boxplot3(2,:)==4) = S_post_wake_boxplot3(3,S_post_wake_boxplot3(2,:)==4)+0.6;
pos = repelem(1:16,1,4);
pos(2:4:end) = pos(2:4:end)+0.2;
pos(3:4:end) = pos(3:4:end)+0.4;
pos(4:4:end) = pos(4:4:end)+0.6;

f7b=figure
hold on
boxplot(S_post_wake_boxplot3(1,:),S_post_wake_boxplot3(3,:),'Symbol','','Positions',pos)
plot((1:width(S_post_wake3))/30/60 + 1,mean(S_post_wake3(id_4g==1 & all_errs3'<=1.2,:),1),'color',c_hi1,'LineWidth',2)
plot((1:width(S_post_wake3))/30/60 + 1,mean(S_post_wake3(id_4g==2 & all_errs3'<=1.2,:),1),'color',c_hi2,'LineWidth',2)
plot((1:width(S_post_wake3))/30/60 + 1,mean(S_post_wake3(id_4g==3 & all_errs3'<=1.2,:),1),'color',c_hi_nlst,'LineWidth',2)
plot((1:width(S_post_wake3))/30/60 + 1,mean(S_post_wake3(id_4g==4 & all_errs3'<=1.2,:),1),'color',c_tem,'LineWidth',2)
h = findobj(gca,'Tag','Box')
for ii=1:4:64; p4=patch(get(h(ii),'XData'),get(h(ii),'YData'),c_tem,'FaceAlpha',.7);end
for ii=2:4:64; p3=patch(get(h(ii),'XData'),get(h(ii),'YData'),c_hi_nlst,'FaceAlpha',.7);end
for ii=3:4:64; p2=patch(get(h(ii),'XData'),get(h(ii),'YData'),c_hi2,'FaceAlpha',.7);end
for ii=4:4:64; p1=patch(get(h(ii),'XData'),get(h(ii),'YData'),c_hi1,'FaceAlpha',.7);end
xticks((1:16)+0.3)
xticklabels(string(1:16))
xlabel('Wake duration (hours)')
ylabel('Theoretical sleep pressure')
lgd=legend([p1 p2 p3 p4],{'IH-LST-S','IH-LST-nS','IH-nLST','HC'},'NumColumns',4,'location','southoutside')
title('Simulated sleep pressure during wakefulness')
grid on

% Stats per hour
pvalues_post_wake3 = [];
for ii=1:16
    val = S_post_wake_boxplot3(1,floor(S_post_wake_boxplot3(3,:))==ii);
    group = S_post_wake_boxplot3(2,floor(S_post_wake_boxplot3(3,:))==ii);
    age = br_list.Age(id_4g>0 & all_errs3'<=1.2);
    sex = br_list.SEXE(id_4g>0 & all_errs3'<=1.2);
    h = lillietest(val(group==1)) + lillietest(val(group==2)) + lillietest(val(group==3)) +lillietest(val(group==4));
    if h>0
        X = [age sex];
        b = regress(val', [ones(size(X,1),1), X]);
        residuals = val' - [ones(size(X,1),1), X]*b;
        [p,tab,stats] = kruskalwallis(residuals,group,'off');
        c = multcompare(stats,'display','off');
        parametric = 0;
    else
        X = [age sex];
        b = regress(val', [ones(size(X,1),1), X]);
        residuals = val' - [ones(size(X,1),1), X]*b;
        [p,tab,stats] = anova1(residuals,group,'off');
        % Adjust
        c=multcompare(stats,'display','off');
        parametric = 1;
    end
    pvalues_post_wake3 = [pvalues_post_wake3,[p;parametric;c(:,6)]];
end


%% FIGURE S1 

mat_theta = {};
id_theta = [];
heure_theta = [];
for ii=find(id_tem_included | id_hi_included)
    tab = spectral_tables{ii};
    hypno = tab.Hypno;
    cycles = tab.Cycles;
    if ismember(ii, id_c4_m2) 
        swa = (tab{:,17}+tab{:,18});
        theta = tab{:,19};
        alpha = tab{:,20};
        beta = tab{:,21} + tab{:,22};
        tot_pwr = sum(tab{:,[17 18 19 20 21 22]},2,'omitnan');
    else
        swa = (tab{:,3}+tab{:,4});
        theta = tab{:,5};
        alpha = tab{:,6};
        beta = tab{:,7} + tab{:,8};
        tot_pwr = sum(tab{:,[3 4 5 6 7 8]},2,'omitnan');
    end
    
    ratio = theta./tot_pwr;
    ratio(cycles>0) = NaN;
    % Remove artefacts : highest 5% of derivative of signal
    dratio = [0;diff(ratio)];    
    artef = dratio>prctile(dratio,95);
    artef = conv(artef,ones(1,3),'same')>=1;
    ratio(artef) = NaN;
    ratio = movmedian(ratio,20,'omitnan');
    
    % Now process during each "long" wake period (>1h30)
    id1 = finddatagroups(hypno,1);
    for jj = 1:2:length(id1)
       if id1(jj+1)-id1(jj)+1 > 60*30
           ratio_cycle = ratio(id1(jj):id1(jj+1));
           id_theta(end+1) = ii;
           mat_theta{end+1} = ratio_cycle;
           heure_theta(end+1) = 23 + (id1(jj))/30/60;
       end
    end
end

%Combine into 1 matrix
mat_theta_aligned = nan(length(mat_theta),max(cellfun(@length,mat_theta)));
mat_theta_lm_aligned = nan(length(mat_theta),max(cellfun(@length,mat_theta)));
for ii=1:length(mat_theta)
   tmp =  mat_theta{ii};  
   mat_theta_aligned(ii,1:length(tmp)) = tmp;
   tmp =  mat_theta_lm{ii};  
   mat_theta_lm_aligned(ii,1:length(tmp)) = tmp;
end

% Now we extract from each wake bout mean relativetheta power over first 
% and last 45 minutes, as well as bout duration

mat_theta2 = [];
for ii=1:height(mat_theta_aligned)
    tmp = mat_theta_aligned(ii,:);
    if sum(~isnan(tmp)) < 90*30
        mat_theta2(ii,:) = [NaN, NaN, NaN];
        continue
    end
    val_start = mean(tmp(1:45*30));
    id = find(~isnan(tmp),1,'last');
    val_end = mean(tmp(id-45*30:id));
    dura = id/30/60;
    mat_theta2(ii,:) = [val_start,val_end,dura];
end

% Stats GLME

mat_theta_g1 = mat_theta2(ismember(id_theta, find(id_4g==1)),:);
timemat1 = mat_theta_g1(:,3);
boutmat1 = repelem([1:height(mat_theta_g1)]',1,1);
sexe1 = repelem(br_list.SEXE(id_theta(ismember(id_theta, find(id_4g==1)))),1,1);
age1 = repelem(br_list.Age(id_theta(ismember(id_theta, find(id_4g==1)))),1,1);
heure1 = heure_theta(ismember(id_theta, find(id_4g==1)));
tbl1 = table(mat_theta_g1(:,2)-mat_theta_g1(:,1),timemat1(:),boutmat1(:),categorical(sexe1(:)),age1(:),heure1(:),'VariableNames',{'Delta','Duration','Bout','Sex','Age','Time'});
mdl1 = fitglme(tbl1,'Delta ~ Duration + Age + Sex + Time + (1|Bout)');

mat_theta_g2 = mat_theta2(ismember(id_theta, find(id_4g==2)),:);
timemat2 = mat_theta_g2(:,3);
boutmat2 = repelem([1:height(mat_theta_g2)]',1,1);
sexe2 = repelem(br_list.SEXE(id_theta(ismember(id_theta, find(id_4g==2)))),1,1);
age2 = repelem(br_list.Age(id_theta(ismember(id_theta, find(id_4g==2)))),1,1);
heure2 = heure_theta(ismember(id_theta, find(id_4g==2)));
tbl2 = table(mat_theta_g2(:,2)-mat_theta_g2(:,1),timemat2(:),boutmat2(:),categorical(sexe2(:)),age2(:),heure2(:),'VariableNames',{'Delta','Duration','Bout','Sex','Age','Time'});
mdl2 = fitglme(tbl2,'Delta ~ Duration + Age + Sex + Time + (1|Bout)');

mat_theta_g3 = mat_theta2(ismember(id_theta, find(id_4g==3)),:);
timemat3 = mat_theta_g3(:,3);
boutmat3 = repelem([1:height(mat_theta_g3)]',1,1);
sexe3 = repelem(br_list.SEXE(id_theta(ismember(id_theta, find(id_4g==3)))),1,1);
age3 = repelem(br_list.Age(id_theta(ismember(id_theta, find(id_4g==3)))),1,1);
heure3 = heure_theta(ismember(id_theta, find(id_4g==3)));
tbl3 = table(mat_theta_g3(:,2)-mat_theta_g3(:,1),timemat3(:),boutmat3(:),categorical(sexe3(:)),age3(:),heure3(:),'VariableNames',{'Delta','Duration','Bout','Sex','Age','Time'});
mdl3 = fitglme(tbl3,'Delta ~ Duration + Age + Sex + Time + (1|Bout)');

mat_theta_g4 = mat_theta2(ismember(id_theta, find(id_4g==4)),:);
timemat4 = mat_theta_g4(:,3);
boutmat4 = repelem([1:height(mat_theta_g4)]',1,1);
sexe4 = repelem(br_list.SEXE(id_theta(ismember(id_theta, find(id_4g==4)))),1,1);
age4 = repelem(br_list.Age(id_theta(ismember(id_theta, find(id_4g==4)))),1,1);
heure4 = heure_theta(ismember(id_theta, find(id_4g==4)));
tbl4 = table(mat_theta_g4(:,2)-mat_theta_g4(:,1),timemat4(:),boutmat4(:),categorical(sexe4(:)),age4(:),heure4(:),'VariableNames',{'Delta','Duration','Bout','Sex','Age','Time'});
mdl4 = fitglme(tbl4,'Delta ~ Duration + Age + Sex + Time + (1|Bout)');

% Now plot

fs1 = figure
hold on
patch([0 15 15 0],[0 0 0.2 0.2], [0.9137    0.9882    0.8627],'EdgeColor','none')
patch([0 15 15 0],[0 0 -0.15 -0.15], [1.0000    0.8588    0.8863],'EdgeColor','none')
s1=scatter(mat_theta2(ismember(id_theta, find(id_4g==1)),3),mat_theta2(ismember(id_theta, find(id_4g==1)),2)-mat_theta2(ismember(id_theta, find(id_4g==1)),1),15,c_hi1,'filled')
s2=scatter(mat_theta2(ismember(id_theta, find(id_4g==2)),3),mat_theta2(ismember(id_theta, find(id_4g==2)),2)-mat_theta2(ismember(id_theta, find(id_4g==2)),1),15,c_hi2,'filled')
s3=scatter(mat_theta2(ismember(id_theta, find(id_4g==3)),3),mat_theta2(ismember(id_theta, find(id_4g==3)),2)-mat_theta2(ismember(id_theta, find(id_4g==3)),1),15,c_hi_nlst,'filled')
s4=scatter(mat_theta2(ismember(id_theta, find(id_4g==4)),3),mat_theta2(ismember(id_theta, find(id_4g==4)),2)-mat_theta2(ismember(id_theta, find(id_4g==4)),1),15,c_tem,'filled')
Xfit = 1.5:0.5:15;
ylin1 =  Xfit.* mdl1.Coefficients.Estimate(2) + mdl1.Coefficients.Estimate(1);
plot(Xfit, ylin1,'color',c_hi1,'LineWidth',3)
ylin2 = Xfit .* mdl2.Coefficients.Estimate(2) + mdl2.Coefficients.Estimate(1);
plot(Xfit, ylin2,'color',c_hi2,'LineWidth',3)
ylin3 = Xfit .* mdl3.Coefficients.Estimate(2) + mdl3.Coefficients.Estimate(1);
plot(Xfit, ylin3,'color',c_hi_nlst,'LineWidth',3)
ylin4 = Xfit .* mdl4.Coefficients.Estimate(2) + mdl4.Coefficients.Estimate(1);
plot(Xfit, ylin4,'color',c_tem,'LineWidth',3)
xlabel('Duration of wake bout (hours)')
ylabel('Difference in relative theta power')
title('Difference in mean theta power in the first and last 45minutes of wake bouts')
xlim([1 15])
legend([s1 s2 s3 s4],{'IH-LSTa','IH-LSTb','IH-nLST','HC','Location','northeast'})

pvals_all = [mdl1.Coefficients.pValue(2);...
    mdl2.Coefficients.pValue(2);...
    mdl3.Coefficients.pValue(2);...
    mdl4.Coefficients.pValue(2)];
pvals_all_adj = mafdr(pvals_all,'BHFDR','true');

