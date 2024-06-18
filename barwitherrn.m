function varargout = barwitherrn(data,groups),

%   Provide data as a matrix (each column will be treated separately), or
%   as grouped data with barwitherr(data, 'grouped'), in which case 'data'
%   is a 2-column matrix with the first column containing the data adn the
%   second column containing group number.
%
%   Original by Martina F. Callaghan; I have changed it quite a bit to now
%   compute means and errors by itself from the data provided, as well as to
%   do tests for significance between groups and from zero for each group.
%
%   This function assumes normality and therefore plots means and standard
%   errors of the mean, and performs parametric significance tests (ttest 
%   for each group to zero and anova between groups). 
%   For a non-parametric equivalent, use <a href="matlab:help barwitherr">barwitherr</a>.
%
%   Raly Todorova

grouped = false;
data = double(data);

if nargin==1,
    if sum(sum(~isnan(data))==0) > 0, %if one of the vectors contains no real values
        varargout = [];
        return
    end
    xOrder = 1:size(data,2);
    values = nanmean(data);
    errors = nansem(data); % - values;
else
    grouped = true;
    if isastring(groups,'groups','grouped','group'),
        groups = data(:,end);
        data = data(:,1:end-1);
    end
    xOrder = unique((groups))';
    n = numel(xOrder);
    values = nan(size(data,2),n);
    errors = nan(size(data,2),n);
    for i=1:n,
        group = xOrder(i);
        values(:,i) = nanmean(data(groups==group,1:end));
        errors(:,i) = nansem(data(groups==group,1:end))';% - values(:,i);
    end
    if size(data,2)>1,
        xOrder = 1:size(data,2);
    end
end

[nRows nCols] = size(values);
if sum(size(xOrder)~=size(values))==0,
    handles.bar = bar(xOrder, values); % standard implementation of bar fn
else
    handles.bar = bar(values); % standard implementation of bar fn
end
hold on
hBar = handles.bar;

if nRows > 1,
    hErrorbar = zeros(1,nCols);
    for col = 1:nCols
        % Extract the x location data needed for the errorbar plots:
%         x = get(get(handles.bar(col),'children'),'xdata');
        x = get(handles.bar(col),'xdata');
        % Use the mean x values to call the standard errorbar fn; the
        % errorbars will now be centred on each bar; these are in ascending
        % order so use xOrder to ensure y values and errors are too:
        hErrorbar(col) = errorbar(mean(x,1), values(xOrder,col), errors(xOrder,col), errors(xOrder, col), '.k');
        set(hErrorbar(col), 'marker', 'none')
    end
else
    x = get(handles.bar,'xdata');
    hErrorbar = errorbar(mean(x,1), values, errors, errors, '.k');
    set(hErrorbar, 'marker', 'none')
end



%% Adding significance indication (stars etc)
if ~grouped || size(data,2)==1,
    if grouped,
        [~,~,stats] = anova1(data(:,1),groups,'off');
        % For each group, check if it's different from zero
        for i=1:unique(groups)';
            thesedata = data(groups==i,1);
            if sum(~isnan(thesedata))>2,
                p = signrank(thesedata);
                if p<0.05,
                    % Put a little star above it
                    plot(xOrder(i), sign(values(i))*(max([abs(values(i)+errors(i)) abs(values(i)-errors(i))])+abs(values(i)/7)), '*', 'markersize', 5, 'MarkerEdgeColor', [0 0 0.5]);
                end
            end
        end
    elseif ~grouped,
        [~,~,stats] = anova1(data,[],'off');
        for i=1:size(data,2),
            thesedata = data(:,i);
            if sum(~isnan(thesedata))>2,
                p = signrank(thesedata);
                if p<0.05, %Put a little star above it to show it's significantly different from zero
                    plot(xOrder(i), sign(values(i))*(max([abs(values(i)+errors(i)) abs(values(i)-errors(i))])+abs(values(i)/7)), '*', 'markersize', 5, 'MarkerEdgeColor', [0 0 0.5]);
                end
            end
        end
    end
    comparison = multcompare(stats,'display', 'off');
    comparison = [comparison(:,1:2) double(comparison(:,3).*comparison(:,5)>0)]; %if the upper and lower bound have the same sign
    comparison2 = multcompare(stats,'display', 'off', 'alpha', 0.01);
    comparison(:,3) = comparison(:,3) + double(comparison2(:,3).*comparison2(:,5)>0); % third column shows the number of stars to be included. 1 for 0.05, 1 more for 0.01, and another one for 0.001       if sum(double(comparison2(:,3).*comparison2(:,5)>0)),
    comparison3 = multcompare(stats,'display', 'off', 'alpha', 0.001);
    comparison(:,3) = comparison(:,3) + double(comparison3(:,3).*comparison3(:,5)>0);
    
    % Now to plot them:
    for i=1:size(comparison,1),
        s = sign(Portion(data(:)>0)./(Portion(data(:)>0 | data(:)<0))-0.5);
        smallnumber = nanmean(values(:))/10; %for display purposes, so things don't overlap
        y1 = s*max(abs([errors(:)+values(:); -errors(:)+values(:)]))+(comparison(i,2)-comparison(i,1))*smallnumber;
        y2 = y1+smallnumber/4;
        x1 = xOrder(comparison(i,1)) - (xOrder(comparison(i,2))-xOrder(comparison(i,1)))*0.2 + 0.4;
        x2 = xOrder(comparison(i,2)) + (xOrder(comparison(i,2))-xOrder(comparison(i,1)))*0.2 - 0.4; %lines between adjecent bars are shorter; the bigger the line, the farther it goes.
        if comparison(i,3)==0,
            %             text(mean([comparison(i,1) comparison(i,2)]), y2, 'n.s.');
        elseif comparison(i,3)==1,
            plot([x1 x2], [y1 y1], 'k');
            plot([x1 x1], [y1 y1-smallnumber/10], 'k');
            plot([x2 x2], [y1 y1-smallnumber/10], 'k');
            text(mean([xOrder(comparison(i,1)) xOrder(comparison(i,2))]), y2, '*');
        elseif comparison(i,3)==2,
            plot([x1 x2], [y1 y1], 'k');
            plot([x1 x1], [y1 y1-smallnumber/10], 'k');
            plot([x2 x2], [y1 y1-smallnumber/10], 'k');
            text(mean([xOrder(comparison(i,1)) xOrder(comparison(i,2))]), y2, '**');
        elseif comparison(i,3)==3,
            plot([x1 x2], [y1 y1], 'k');
            plot([x1 x1], [y1 y1-smallnumber/10], 'k');
            plot([x2 x2], [y1 y1-smallnumber/10], 'k');
            text(mean([xOrder(comparison(i,1)) xOrder(comparison(i,2))]), y2, '***');
        end
    end
    
else
    x = get(handles.bar,'children');
    for i=1:n,
        xs(:,i) = nanmean(get(x{i},'xdata'));
    end
    for j=1:length(xOrder),
        [~,~,stats] = anova1(data(:,j), groups,'off');
        comparison = multcompare(stats,'display', 'off');
        comparison = [comparison(:,1:2) double(comparison(:,3).*comparison(:,5)>0)]; %if the upper and lower bound have the same sign
        comparison2 = multcompare(stats,'display', 'off', 'alpha', 0.01);
        comparison(:,3) = comparison(:,3) + double(comparison2(:,3).*comparison2(:,5)>0); % third column shows the number of stars to be included. 1 for 0.05, 1 more for 0.01, and another one for 0.001       if sum(double(comparison2(:,3).*comparison2(:,5)>0)),
        comparison3 = multcompare(stats,'display', 'off', 'alpha', 0.001);
        comparison(:,3) = comparison(:,3) + double(comparison3(:,3).*comparison3(:,5)>0);
        sigFor2Groups(j,1) = comparison(1,3);
        for i=1:size(comparison,1),
            s = sign(Portion(data(:,j)>0)./(Portion(data(:,j)>0 | data(:,j)<0))-0.5);
            smallnumber = nanmean(values(:))/10; %for display purposes, so things don't overlap
            y1 = s*max(abs([errors(j,:)'+values(j,:)'; -errors(j,:)'+values(j,:)']))+(comparison(i,2)-comparison(i,1))*smallnumber;
            y2 = y1+smallnumber/3;
            x1 = xs(j,comparison(i,1));
            x2 = xs(j,comparison(i,2));
            if comparison(i,3)==1,
                plot([x1 x2], [y1 y1], 'k');
                plot([x1 x1], [y1 y1-smallnumber/10], 'k');
                plot([x2 x2], [y1 y1-smallnumber/10], 'k');
                text(mean([x1 x2]), y2, '*','HorizontalAlignment','center');
            elseif comparison(i,3)==2,
                plot([x1 x2], [y1 y1], 'k');
                plot([x1 x1], [y1 y1-smallnumber/10], 'k');
                plot([x2 x2], [y1 y1-smallnumber/10], 'k');
                text(mean([x1 x2]), y2, '**','HorizontalAlignment','center');
            elseif comparison(i,3)==3,
                plot([x1 x2], [y1 y1], 'k');
                plot([x1 x1], [y1 y1-smallnumber/10], 'k');
                plot([x2 x2], [y1 y1-smallnumber/10], 'k');
                text(mean([x1 x2]), y2, '***','HorizontalAlignment','center');
            end
        end
    end
    comparison = sigFor2Groups;
end





hold off

switch nargout
    case 1
        varargout{1} = hBar;
    case 2
        varargout{1} = hBar;
        varargout{2} = hErrorbar;
    case 3
        varargout{1} = hBar;
        varargout{2} = hErrorbar;
        varargout{3} = comparison;
end


