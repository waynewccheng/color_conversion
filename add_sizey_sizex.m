%% add sizex and sizey into reflectance.mat and rename to transmittance
% for WSI truthing project
% change name from reflectance to transmittance
% 11-24-2018 Thanksgiving
% 4-9-2018

function add_sizey_sizex (folder_name)

    err = load([folder_name '\reflectance.mat'],'sizey');

    already_added = any(structfun(@isempty,err))
    
    if ~already_added
       load([folder_name '\vimarray.mat'],'vimarray');
       [wl, sizey, sizex] = size(vimarray)
    end
    
    save([folder_name '\reflectance.mat'],'sizey','sizex','-append')

end
