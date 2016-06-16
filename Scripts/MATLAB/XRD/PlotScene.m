function PlotScene(Scene,N,M)

if ~exist('N')
  N=1;  % First scene
end
if ~exist('M')
  M=1;  % Type of plot
end
Nmax=length(Scene);

Screen=get(0,'Screensize');
figure('color','w','position',[Screen(3)*0.2 Screen(4)*0.2 Screen(3)*0.6 Screen(4)*0.6])
cmap=colormap(lines(Nmax));

% Initial Plot
plotting(N,M)


popup = uicontrol('Style', 'popup',...
  'String', {'Color','Black & White'},...
  'Position', [70 340 100 50],...
  'Callback', @setmap,'Value',M);

sld = uicontrol('Style', 'slider',...
  'Min',1,'Max',Nmax,'Value',N,...
  'Position', [180 70 500 40],...
  'Callback', @scene_select);


function scene_select(source,callbackdata)
  N=round(sld.Value);
  sld.Value=N;
  plotting(N,M)
end

function setmap(source,callbackdata)
  val = source.Value;
  maps = source.String;
  type = maps{val};
  if strcmpi(type,'Black & White')
    M=2;
  else 
    M=1;
  end
  plotting(N,M)
end

function plotting(N,M)
  subplot(1,2,1)
  if M==1;
    imshow(Scene(N).Image)
  else
    imshow(Scene(N).BW); hold on
    
    % Plot edges
    for k = 1:length(Scene(N).Edges)
      boundary = Scene(N).Edges{k};
      plot(boundary(:,2), boundary(:,1), 'color',cmap(N,:), 'LineWidth', 2)
    end
      
    % Plot centroids
    Centers = Scene(N).Centers;
    text(Centers(:,1),Centers(:,2),num2str(N),'color',cmap(N,:))
    
  end
  title(Scene(N).Filename,'interpret','none')
    
  subplot(1,2,2)
  set(gca, 'YDir', 'reverse');
  [a,b,~]=size(Scene(N).Image);
  axis([0 a 0 b]);
end



subplot(1,2,2); hold on
for i=1:length(Scene)
  % Draw edges
  for k = 1:length(Scene(i).Edges)
    boundary = Scene(i).Edges{k};
    plot(boundary(:,2), boundary(:,1), 'color',cmap(i,:), 'LineWidth', 2)
  end
  
  % Plot centroids
  Centers = Scene(i).Centers;
  text(Centers(:,1),Centers(:,2),num2str(i),'color',cmap(i,:))
  
end
set(gca, 'YDir', 'reverse');
[a,b,~]=size(Scene(N).Image);
axis([0 a 0 b]);
axis square

end