figure;subplot(3,1,1);
vals=[optomeanMat{3,2}(2,1) -1*(optomeanMat{3,2}(2,1)-1) ...
    -1*(optomeanMat{3,3}(2,1)-1) (optomeanMat{3,3}(2,1))];
vals=vals*100;
b=bar(vals);
text(b.XEndPoints, b.YEndPoints, string(vals), ...
    VerticalAlignment='bottom', HorizontalAlignment='center')
hold on; title('Reinforcement');
xticklabels({'Hit','Miss','CR','FA'});

subplot(3,1,2);
vals=[optomeanMat{3,4}(2,1) -1*(optomeanMat{3,4}(2,1)-1) ...
    -1*(optomeanMat{3,5}(2,1)-1) optomeanMat{3,5}(2,1)];
vals=vals*100;
b=bar(vals);text(b.XEndPoints, b.YEndPoints, string(vals), ...
    VerticalAlignment='bottom', HorizontalAlignment='center')
hold on; title('Fifth Dead');
xticklabels({'Hit','Miss','CR','FA'});

subplot(3,1,3); clear vals
vals=[optomeanMat{3,6}(2,1) -1*(optomeanMat{3,6}(2,1)-1) ...
    -1*(optomeanMat{3,7}(2,1)-1) optomeanMat{3,7}(2,1)];
vals=vals*100;
b=bar(vals);
text(b.XEndPoints, b.YEndPoints, string(vals), ...
    VerticalAlignment='bottom', HorizontalAlignment='center')
hold on; title('Second Dead');
xticklabels({'Hit','Miss','CR','FA'});
ylim([0 100]);