function selective_classifier(varargin)
% Function that classifies data based on user inputs
%
% Inputs:
% ========================================================================
% stim_data        : (2xN) 2 columns of stim hist data to serve as (x,y) coordinates for classification points
% animal_class     : (1xN cell) cell vector containing selectivity classsification text "Above 90" or "Below 90"
% animal_names     : (1xN cell) cell vector containing the name of each animal associated with with stim data and dcassification text
% fig_title        : (1xN cell) figure title
% unit_type        : (char) name of unit type to display and scale data to
% arg1             : (1xN num) vector of data to assist in data scaling
% split_by_animal  : (bool) classify data based on each animal
% show_shared_chans: (bool) show text above points showing the number of
%                           shared channels responding to stim
% shared_chans     : (1xN cell) cell vector of characters showing the number of channels activated corresponding to
%                         each stim  pair
% 
% Example: selective_classifier(data, selectivity_labels, animal_names, fig_title, 'threshold', minthresh_stim, false, true, chan_txt);

% If not splitting by animal, set the "classifier_type" to desired classifier
classifier_type = 'LDA';


data = varargin{1};
stim_data = data(:,[1,3]);
% separate and sort channel responses from cuff pairs
[~,idx] = sort(data(:,1)); % sort just the first 2 columns
data_1_2 = data(idx,1:2);
[~,idx] = sort(data(:,4)); % sort just the last 2 columns
data_3_4 = data(idx,3:4);

animal_class = varargin{2};
animal_names = varargin{3};

if length(varargin) > 3
    fig_title = varargin{4};
else
    fig_title = '';
end

if length(varargin) > 4
    unit_type = varargin{5};
    arg1 = varargin{6};
    switch unit_type
        case 'charge'
            stim_data = (stim_data.*arg1)/1000000;
            units = '(C)';
            x_label = ['Cuff 1-2 Stim', units];
            y_label = ['Cuff 3-4 Stim', units];
        case 'amp'
            min_1_2 = 0;
            min_3_4 = 0;
            units = '(\muA)';
            x_label = ['Cuff 1-2 Stim', units];
            y_label = ['Cuff 3-4 Stim', units];
        case 'threshold'
            % arg1 must be a list of min values to divide data by, since each
            % animal and pw has its own min thresh value
            stim_data = stim_data./arg1;
            min_1_2 = min(stim_data(:,1));
            min_3_4 = min(stim_data(:,2));
            
            x_label = 'Cuff 1-2 Threshold Multiple';
            y_label = 'Cuff 3-4 Threshold Multiple';
    end
else
    unit_type = 'amp';
    units = '(\muA)';
    x_label = ['Cuff 1-2 Stim', units];
    y_label = ['Cuff 3-4 Stim', units];
end

if length(varargin)>6
    split_by_animal = varargin{7};
else
    split_by_animal = true;
end

if length(varargin)>7 && varargin{8}==true
    show_shared_chans = true;
    shared_chans = varargin{9};
else
    show_shared_chans = false;
    shared_chans = '';
end



figure(8);

max_1_2 = max(stim_data(:,1));
max_3_4 = max(stim_data(:,2));
%if max_3_4 > max_1_2
%    max_1_2 = max_3_4;
%else
%    max_3_4 = max_1_2;
%end

% plot data points (correctly and incorrectly classified)
if split_by_animal==true
    NUM_K = numel(unique(animal_names));
    unique_animal_names = cell2mat(unique(animal_names));
    %color_list = {'k','m','c','r','g','b'}; % simpler list
    color_list = {
        [0, 0, 1],...
        [0, 0.5, 0],...
        [1,0,0],...
        [0, 0.75, 0.75],...
        [0.75, 0, 0.75],...
        [0.8500, 0.3250, 0.0980],...
        };
    
    a = 0.1;
    b = 0.5;
    
    if size(unique_animal_names,1) == 1
        subplot(4,4,[2:4,6:8,10:12])
    end
    
    % Better to split by selective points before animals
    temp_idx = find(strcmp(animal_class, 'Above 90'));
    temp_data = stim_data(temp_idx,:);
    temp_animal_names = animal_names(temp_idx);
    unique_stim_vals = [];
    if show_shared_chans==true
        for k=1:length(temp_idx)
            temp_shared_chans{k,1} = shared_chans{temp_idx(k),1};
        end
    end
    
    %Plot data points according to associated animal
    for i=1:size(unique_animal_names,1)
        temp_animal_names_idx = find(strcmp(temp_animal_names, unique_animal_names(i,:)));
        new_data_selec = temp_data(temp_animal_names_idx,:);
        unique_stim_vals = [unique_stim_vals;new_data_selec];
        s = scatter(new_data_selec(:,1), new_data_selec(:,2), 'LineWidth', 0.2, 'MarkerEdgeColor','k','MarkerFaceColor', color_list{i},'Marker','o','SizeData',50 + i*5);
        s.MarkerFaceAlpha = (b-a).*rand(1) + a;
        hold on;
        if show_shared_chans==true
            for j=1:size(new_data_selec,1)
                text(new_data_selec(j,1), new_data_selec(j,2),temp_shared_chans{j,1},'HorizontalAlignment','right','VerticalAlignment','top','FontSize',10);
            end
        end
    end
    
    % Better to split by selective points before animals
    temp_idx = find(strcmp(animal_class, 'Below 90'));
    temp_data = stim_data(temp_idx,:);
    temp_animal_names = animal_names(temp_idx);
    if show_shared_chans==true
        for k=1:length(temp_idx)
            temp_shared_chans{k,1} = shared_chans{temp_idx(k),1};
        end
    end
    
    %Plot data points according to associated animal
    for i=1:size(unique_animal_names,1)
        temp_animal_names_idx = find(strcmp(temp_animal_names, unique_animal_names(i,:)));
        new_data_non_selec = temp_data(temp_animal_names_idx,:);
        unique_stim_vals = [unique_stim_vals; new_data_non_selec];
        s = scatter(new_data_non_selec(:,1), new_data_non_selec(:,2), 'LineWidth', 0.2, 'MarkerEdgeColor','k','MarkerFaceColor', color_list{i}, 'Marker','s','SizeData',25);
        s.MarkerFaceAlpha = (b-a).*rand(1) + a;
        %if show_shared_chans==true
        %    for j=1:size(new_data,1)
        %        text(new_data(j,1), new_data(j,2),temp_shared_chans{j,1},'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',10);
        %    end
        %end
    end
    
    chan_resp_1_2 = [0,0; unique(data_1_2,'rows')];
    chan_resp_3_4 = [0,0; unique(data_3_4,'rows')];
    
    
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    %xlabel(x_label);
    %ylabel(y_label);
    legend(unique_animal_names)
    xlim([min_1_2-min_1_2/5 max_1_2 + 0.5]);
    ylim([min_3_4-min_3_4/5 max_3_4 + 0.5]);
    
    if strcmp(unit_type,'threshold')==1
        chan_resp_1_2(:,1) = chan_resp_1_2(:,1)./arg1(1,1);
        chan_resp_3_4(:,1) = chan_resp_3_4(:,1)./arg1(1,2);
    end
    
    
    % Works only if there is one amimal to display channel data, otherwise
    % it gets over saturated
    if  size(unique_animal_names,1) == 1
        %x-axis plot shows number of channels responding from cuff3-4 stim
        subplot(4,4,(14:16))
        plot(chan_resp_1_2(:,1),chan_resp_1_2(:,2));
        xlabel(x_label)
        ylabel('# Chans');
        xlim([min_1_2-min_1_2/5 max_1_2 + 0.5])
        ylim([0 32])
        %y-axis plot shows number of channels responding from cuff1-4 stim
        subplot(4,4,([1,5,9]))
        plot(chan_resp_3_4(:,1),chan_resp_3_4(:,2));
        set(gca,'XAxisLocation','Top');
        xlabel(x_label)
        camroll(90)
        ylabel('# Chans');
        ylim([0 32])
        xlim([min_3_4-min_3_4/5 max_3_4 + 0.5])
        
        sgtitle(fig_title);
    else
        title(fig_title)
    end
    hold off;
    drawnow;
else
    NUM_K = ['Above 90','Below 90'];
    for k = 1:size(animal_class,1)
        if animal_class{k,1}=='Above 90'
            num_animal_class(k,1) = 1;
        else
            num_animal_class(k,1) = 0;
        end
    end
    
    % # plot grid classification color-coded
    % hold on
    % image(X, Y, reshape(grp2idx(C),npoints,npoints))
    % axis xy, colormap(clrLite)
    numInst = size(stim_data,1);             %# number of instances
    
    % Create mesh grid data points for classifier shading
    npoints = 50;
    mn = min(stim_data);    mx = max(stim_data);
    [X,Y] = meshgrid( linspace(mn(1),mx(1),npoints) , linspace(mn(2),mx(2),npoints) );
    X = X(:); Y = Y(:);
    
    switch classifier_type
        
        case 'fisher_LDA'
            
            
            %[w,t,fp]=fisher_training(stim_data,num_animal_class);
            %[l,precision,recall,accuracy,F1]=fisher_testing(stim_data,w,t,num_animal_class);
            %classError = 1-accuracy;
            %xx=0:0.1:max([max_1_2,max_3_4]);
            %yy=-w(1)/w(2)*xx+t/w(2);
            %plot(xx,yy,'-.k','LineWidth',1);
            
        case 'LDA'
            % ==================  Linear Discriminant Analysis  ===================
            mdl = fitcdiscr(stim_data, num_animal_class);
            classError = loss(mdl, stim_data, num_animal_class);
            lx = 0:0.1:(max([max_1_2,max_3_4])+0.5);
            %lx = get(gca, 'Xlim');
            ly = -(mdl.Coeffs(1, 2).Const + mdl.Coeffs(1, 2).Linear(1) .* lx) / mdl.Coeffs(1, 2).Linear(2);
            %plot(lx, ly, '-.k', 'DisplayName', 'LDA')
            %[C,err,P,logp,coeff] = classify([X Y], stim_data, animal_class, 'linear');
            % # find incorrectly classified training data
            %[CPred,err] = classify(stim_data(:,1), stim_data(:,2), animal_class, 'linear');
            %bad = ~strcmp(CPred,animal_class);
            %if ~exist('bad','var')
            %    bad = logic(zeros(length(animal_class)));
            %end
            %classError = sum(bad)/numInst;
            
            % draw decision boundaries between pairs of clusters
            %K = coeff(1,2).const;
            %L = [coeff(1,2).linear(1), coeff(1,2).linear(2)];
            %f = sprintf('0 = %g + %g*x + %g*y', K,L(1),L(2));
            %h2 = ezplot(f, [0 max_1_2 0 max_3_4]);
            %set(h2, 'LineStyle', '-.','Color','k', 'LineWidth',1,{'DisplayName'},{'LDA Boundary'})
            
        case 'QDA'
            % ================= Quadratic Discriminant Analysis ================
            [C,err,P,logp,coeff] = classify([X Y], stim_data, animal_class, 'quadratic');
            % # find incorrectly classified training data
            [CPred,err] = classify(stim_data(:,1), stim_data(:,2), animal_class, 'quadratic');
            bad = ~strcmp(CPred,animal_class);
            if ~exist('bad','var')
                bad = logic(zeros(length(animal_class)));
            end
            classError = sum(bad)/numInst;
            
            K = coeff(1,2).const;
            L = coeff(1,2).linear;
            Q = coeff(1,2).quadratic;
            %Function to compute K + L*v + v'*Q*v for multiple vectors
            %v=[x;y]. Accepts x and y as scalars or column vectors.
            f = @(x,y) K + L(1)*x + L(2)*y + Q(1,1)*x.*x + (Q(1,2)+Q(2,1))*x.*y + Q(2,2)*y.*y;
            h2 = fimplicit(f,[0 2000 0 2000]);
            set(h2,'Color','m','LineWidth',2,'DisplayName','Decision Boundary')
            
        case 'naive'
            % =====================  Region Classifier  =====================
            % Naive Bayes
            cp = cvpartition(animal_class,'KFold',10);
            nbGau = fitcnb(stim_data,animal_class);
            nbGauResubErr = resubLoss(nbGau);
            nbGauCV = crossval(nbGau, 'CVPartition',cp);
            nbGauCVErr = kfoldLoss(nbGauCV);
            labels = predict(nbGau, [X Y]);
            gscatter(X,Y,labels,'rb','so')
            classError = 1-nbGauResubErr;
        case 'tree'
            % =======================  Decision Tree  =========================
            t = fitctree(stim_data(:,1:2), animal_class,'PredictorNames',{'SL' 'SW' });
            [grpname,node] = predict(t,[X Y]);
            gscatter(X,Y,grpname,'grb','sod')
            
        case 'kNN'
            % ==========================  kNN  ============================
            chiSqrDist = @(x,Z,wt)sqrt((bsxfun(@minus,x,Z).^2)*wt);
            fprintf("\n k = ")
            for k=1:40
                fprintf("%d, ",k)
                %k = input('Choose value for k: ');
                w = [0.75; 0.25];
                KNNMdl = fitcknn(stim_data,animal_class,'Distance',@(x,Z)chiSqrDist(x,Z,w),...
                    'NumNeighbors',k,'Standardize',1);
                rng(1); % For reproducibility
                CVKNNMdl = crossval(KNNMdl);
                classError = kfoldLoss(CVKNNMdl);
                err_list(k) = classError;
                [grpname,node] = predict(KNNMdl,[X Y]);
                %                gscatter(X,Y,grpname,'grb','sod');
                %                hold on
                %                gscatter(stim_data(:,1), stim_data(:,2), animal_class, clrDark, '.', 20, 'on');
                %                plot(NaN,NaN,'.k');
                %                legend('Above 90%','Below 90%',sprintf('accuracy = %.2f%%', 100*(1-classError)),'Location','northeast')
            end
            [~, k_final] = max(1-err_list);
            fprintf("Best k = %d", k_final)
            KNNMdl = fitcknn(stim_data,animal_class,'Distance',@(x,Z)chiSqrDist(x,Z,w),...
                'NumNeighbors',k_final,'Standardize',1);
            rng(1); % For reproducibility
            CVKNNMdl = crossval(KNNMdl);
            classError = kfoldLoss(CVKNNMdl);
            
            [grpname,node] = predict(KNNMdl,[X Y]);
            gscatter(X,Y,grpname,'grb','sod');
    end
    
    %hold on
    % Better to split by selective points before animals
    temp_idx = find(strcmp(animal_class, 'Above 90'));
    temp_above_data = stim_data(temp_idx,:);
    scatter(temp_above_data(:,1), temp_above_data(:,2), 20, 'MarkerEdgeColor', [0 0 0.8], 'MarkerFaceColor', [0 0 0.8]);
    hold on
    temp_idx = find(strcmp(animal_class, 'Below 90'));
    temp_below_data = stim_data(temp_idx,:);
    scatter(temp_below_data(:,1), temp_below_data(:,2), 20, 'MarkerEdgeColor', [0.8 0.8 0], 'MarkerFaceColor', [0.8 0.8 0]);
    % Plot data points and error
    
    %clrDark = [0 0 0.8 ; 0.8 0.8 0 ; 0.6 0.1 0.6 ; 0 0 0.8 ; 0.5 0 0.8 ; 0 0.5 0.8];
    if ~exist('classError')
        classError = 0;
    end
    legend('Above 90%','Below 90%');
    %if length(NUM_K) > 1
        %plot(NaN,NaN,'.k');
    %    legend(sprintf('Accuracy = %.2f%%', 100*(1-classError)),'Selective','Non-Selective','Location','northeast')
        %legend('Below 90%','Above 90%',sprintf('accuracy = %.2f%%', 100*(1-sum(bad)/numInst)),'Location','northeast')
    %end
    hold off;
    xlim([0 max_1_2 + 0.5]);
    ylim([0 max_3_4 + 0.5]);
    title([fig_title, ' Classifier: ', classifier_type])
    xlabel(x_label);
    ylabel(y_label);
    drawnow;
end


%title( sprintf('accuracy = %.2f%%', 100*(1-sum(bad)/numInst)) )
%xlim([xmin-(xmin/2) xmax+(xmax/20)]);
%ylim([ymin-(ymin/2) ymax+(ymax/20)]);



end


function num_list = shorten_data(data, cuff_pair)
num_list = [];
switch cuff_pair
    case '1-2'
        for i=1:(length(data)-1)
            if data(i)~=data(i+1)
                num_list = [num_list, data(i)];
            end
        end
    case '3-4'
        for i=1:(length(data)-1)
            if ~ismember(data(i), num_list)
                
            end
        end
end

end

