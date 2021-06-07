%% Optional inputs the default values are shown if the input is
clear all
f = figure;
ax = gca;
opts.LineStyle='none';
resultsPath = '../results/';

%% Generate a DOE
DOE = createDOE();

%% Example loop to repeat for all
for i = 1:length(DOE)
data = DOE(i);
data.ID = i;

data.name = sprintf('%01.0f_%s.png', [data.ID, data.equation]);
data.method = "implicit";
data.tform = transformObject(data); %Create transform object
%% Run the script
data = main(data);


%% handle the plotting
if data.metrics.errorFlag ==0
    cla(ax,'reset');
    axis(ax,"equal","vis3d");  view(ax,3);
    plotFigure(ax,'tri-mesh',opts,data);
    caxis(ax,[0 180]);
    
    %% handle exporting figures
    fid = strcat(resultsPath,data.name);
    saveas(f,fid);
    
    %% handle outputs for DOE analysis
    DOEout(i) = data.metrics;
end
end

%% Generate examples for alternative style plots
    cla(ax,'reset');
    axis(ax,"equal","vis3d");  view(ax,3);
    plotFigure(ax,'pole-fig',opts,data);
    fid = strcat(resultsPath,'Pole_FigureExample.png');
    saveas(f,fid);

    cla(ax,'reset');
    axis(ax,"equal","vis3d");  view(ax,3);
    plotFigure(ax,'field-viewer',opts,data);
    fid = strcat(resultsPath,'Example_FieldViewer.png');
    saveas(f,fid);

    cla(ax,'reset');
    plotFigure(ax,'histogram',opts,data);
    fid = strcat(resultsPath,'Example_Histogram.png');
    saveas(f,fid);

    cla(ax,'reset');
    plotFigure(ax,'histogram2',opts,data);
    fid = strcat(resultsPath,'Example_Histogram2.png');
    saveas(f,fid);

    cla(ax,'reset');
    data.p1 = 'surfaceArea'; data.p2 = 'volume';
    plotFigure(ax,'metric-history',opts,data);
    fid = strcat(resultsPath,'Example_MetricHistory.png');
    saveas(f,fid);
