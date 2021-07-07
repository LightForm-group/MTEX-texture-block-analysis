function [  ] = IPF_map_plot(phase, ebsd, outputFileName, visible);
% A function for plotting alpha (hexagonal close-packed) or beta (body-centred-cubic) pole figures using MTEX
  
  if strcmp(visible, 'on');  
    set(0,'DefaultFigureVisible','on'); 
    set(groot,'DefaultFigureVisible','on');
  elseif strcmp(visible, 'off');
    set(0,'DefaultFigureVisible','off');
    set(groot,'DefaultFigureVisible','off');
  else
    disp ('Visibility of pole figures not set as on or off.');
    return;
  end

  if strcmp(phase, 'alpha');
    IPF_map = figure();
    plot(ebsd,ebsd.bc);
    colormap gray
    mtexColorbar;
    hold on
    oM = ipfHSVKey(ebsd('Ti-Hex'))
    oM.inversePoleFigureDirection = yvector;
    color = oM.orientation2color(ebsd('Ti-Hex').orientations);
    plot(ebsd('Ti-Hex'),color);
    hold off
    saveas (IPF_map, outputFileName, 'png');
    close(IPF_map);
    
  elseif strcmp(phase, 'beta');
    IPF_map = figure();
    plot(ebsd);
    hold on
    oM = ipfHSVKey(ebsd('Titanium Cubic'))
    oM.inversePoleFigureDirection = yvector;
    color = oM.orientation2color(ebsd('Titanium Cubic').orientations);
    plot(ebsd('Titanium Cubic'),color);
    hold off
    saveas (IPF_map, outputFileName, 'png');
    close(IPF_map);
    
  else 
    disp ('Phase not recognised for plotting IPF map.');
    return;
  end
  
end