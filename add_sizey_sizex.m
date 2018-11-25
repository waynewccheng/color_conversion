% 11-24-2018
% 4-9-2018
%% add sizex and sizey into the reflectance file for WSI truthing project

function add_sizey_sizex (folder_name)

    err = load([folder_name '\reflectance.mat'],'sizey');

    already_added = any(structfun(@isempty,err))
    
    if ~already_added
       load([folder_name '\vimarray.mat'],'vimarray');
       [wl, sizey, sizex] = size(vimarray)
    end
    
    save([folder_name '\reflectance.mat'],'sizey','sizex','-append')
end
