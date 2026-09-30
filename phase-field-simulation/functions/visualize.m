% Plot temperature/shape
function visualize(T,phi,r,z,Lr,Lz)
    
    global nr %#ok<*GVMIS> 
    
    T = T(1:nr,:);
    Tplt = [flipud(T);T]';
    cla; hold on
   
    cmap = flipud(cmocean('ice'));
    cmap = cmap(1:225,:);
    colormap(cmap);
    
    imagesc('xdata',[-r(nr) r(nr)],'ydata',[0 Lz],'cdata',Tplt,[0 1])
    axis equal tight
    axis([-r(nr) r(nr) 0 Lz])
    box on
    xticks([]); yticks([]);
    set(gca,'xticklabels',[],'yticklabels',[]) 
    
    hb = colorbar;
    hb.Ticks = [0 1]; hb.TickLabelInterpreter = 'latex'; hb.FontSize = 16; 
    hb.Title.String = '$T/T_\infty$'; hb.Title.Interpreter = 'latex';

    fig = gcf;
    fig.Units = 'inches';
    fig.PaperSize = [5 4.5];
    fig.PaperPosition = [0 0 fig.PaperSize];
    fig.Position = [5 5 0 0] + fig.PaperPosition;

end