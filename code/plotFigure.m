function plotFigure(ax,figureType,opts,data)
%updateCellFigs takes input data and plots to axis ax
% Inputs:
%   ax          - handle of axis to plot to
%   figureType  - type of figure to plot
%   (tri-mesh,histogram,pole-fig,histogram2, 
%   opts        - options for plotting
%   data        - structured dataset

opts.EdgeAlpha = 0.1;
opts.LineWidth = 0.5;
opts.FaceColor = 'flat';
n = 100;

% Handle input if p1 or p2 are not defined
if ~isfield(data,'p1')
    data.p1 = "inclination angle (^odeg)";
end
if ~isfield(data,'p2')
    data.p2 = "mean curvature (mm^-^1)";
end

prop1 = extractBefore(data.p1,' ');
prop2 = extractBefore(data.p2,' ');


switch figureType
    case 'tri-mesh'
        % Axis view options, labels & colormap options
        if ~isempty(data.FV.faces)
            patch(ax,'Faces',data.FV.faces,'Vertices',data.FV.vertices,'FaceVertexCData',data.FV.(prop1),opts);
        end
        if ~isempty(data.FVcap.faces)
            hold(ax,"on"); patch(ax,'Faces',data.FVcap.faces,'Vertices',data.FVcap.vertices,'FaceColor',[0.25,0.25,0.25],...
                'FaceAlpha',0.9,'LineStyle','none');
            hold(ax,"off");
        end
        xlabel(ax,"X (mm)"); ylabel(ax,"Y (mm)"); zlabel(ax,"Z - Build Direction (mm)");
        c=colorbar; c.Label.String = data.p1; colormap(ax,"jet"); caxis(ax,'auto');
        
    case 'histogram'
        hold(ax,'on');
        [x, xp, xc] = areaWeightedHistogram(data.FV.(prop1),data.FV.Area,n);
        bar(ax,x,xp); grid(ax,"on"); ylabel(ax,"Probability Density");
        yyaxis(ax,'right'); plot(ax,x,xc,'-r');
        xlabel(ax,data.p1); ylabel(ax,"Cumulative Distribution");
        
    case 'pole-fig'
        [Vs,Fs,Cs] = poleFigure([data.FV.Nx,data.FV.Ny,data.FV.Nz],data.FV.Area,3);
        patch(ax,'Faces',Fs,'Vertices',Vs,'CData',Cs,opts);
        xlabel(ax,"X"); ylabel(ax,"Y"); zlabel(ax,"Z");
        c=colorbar; c.Label.String = "Probability Density"; colormap(ax,hot);
        
    case 'histogram2'
        histogram2(ax,data.FV.(prop1),data.FV.(prop2),n,...
            "DisplayStyle","tile","ShowEmptyBins","on","Normalization","pdf");
        xlabel(ax,data.p1); ylabel(ax,data.p2); 
        c=colorbar; c.Label.String = "Probability Density"; colormap(ax,hot);
    
    case 'field-viewer'
        % Sample Points
        opts.FaceAlpha = 0.5;
        patch(ax,'Faces',data.FV.faces,'Vertices',data.FV.vertices,'FaceVertexCData',[0 0 0.25],...
            'FaceColor','flat','FaceAlpha',0.5,'LineStyle','none'); hold(ax,'on');
        [xq, yq, zq] = meshgrid(data.field.xq,data.field.yq,40/100*data.bulkSize(3));
        xyslice = interp3(data.field.X,data.field.Y,data.field.Z,data.field.(prop1),xq,yq,zq);
        surf(ax,xq,yq,zq,xyslice,'LineStyle','none');
        [S,L] = bounds(data.field.(prop1),'all','omitnan');
        c=colorbar; c.Label.String = data.p1; colormap(ax,jet);
        caxis(ax,'manual'); caxis(ax,[S L]);
        
    case 'metric-history'
        try
            x = data.metrics.(prop1);
            y = data.metrics.(prop2);
        catch
            x = data.metrics.volume;
            y = data.metrics.surfaceArea;
        end
        scatter(ax,x,y,'red',"filled");
        xlabel(ax,data.p1);ylabel(ax,data.p2)
end
end
