function [out] = ginput3d(N)
% GINPUT3D allows you to ginput in xy, yz and zx plane
% of 3D plot, ginput can do that only in xy plane
% In xy plane value of z is set to 0
% In xz plane value of y is set to 0
% In yz plane value of x is set to 0
%% IMPORTANT:
%% work only in xy, yz and xz plane


%Examples:

%[x,y,z] = ginput3d(5);


out=zeros(N,3);
%
if nargin==0
       
error('There are not enought input arguments');
    
elseif nargin>1

error('Too many input arguments');     
        
elseif nargin==1           
     

 % Suspend figure functions
    state = uisuspend(gcf);
    
    toolbar = findobj(allchild(gcf),'flat','Type','uitoolbar');
    if ~isempty(toolbar)
        ptButtons = [uigettool(toolbar,'Plottools.PlottoolsOff'), ...
            uigettool(toolbar,'Plottools.PlottoolsOn')];
        ptState = get (ptButtons,'Enable');
        set (ptButtons,'Enable','off');
    end
 % End of Suspend figure functions
 
 
set(gcf,'Pointer','fullcross')  



    for i=1:N
        
        if ~waitforbuttonpress
        
            XYZ_inic= get(gca,'currentpoint');
            XYZ=~(XYZ_inic(1,:)-XYZ_inic(2,:)).*XYZ_inic(1,:);

            out(i,:)=XYZ;
        
        end

    end

set(gcf,'Pointer','arrow')
refresh 

 % Enable figure functions
       uirestore(state);
    if ~isempty(toolbar) && ~isempty(ptButtons)
        set (ptButtons(1),'Enable',ptState{1});
        set (ptButtons(2),'Enable',ptState{2});
    end
     % End of figure functions



end




