% make json file
subj = input('Enter subject ID (e.g. S0014): ', 's');
user_id = input('Enter user ID (e.g. jsh3653): ', 's');
imaging_dir = sprintf('/Users/%s/Library/CloudStorage/Box-Box/ECoG_imaging/%s_complete/elecs/', user_id, subj);
load(sprintf('%s/TDT_elecs_all.mat', imaging_dir));
json_file = sprintf('%s/%s_electrodes.json', imaging_dir, subj);

fid = fopen(json_file, 'w');

fprintf(fid,'[');
for i=1:size(elecmatrix,1)
    anat = anatomy(i,4);
    anat = anat{1};
    shortname = anatomy(i,1);
    shortname = shortname{1};
    numstart = regexp(shortname, '\d', 'once');
    devname = shortname(1:numstart-1);
    if elecmatrix(i,1) < 0
        hem = 'L';
    else
        hem = 'R';
    end
    if i==size(elecmatrix,1)
        json = fprintf(fid,'{"subid":"%s", "elecid":"%s", "ElecType":"D", "Hem":"%s", "lepto":[%.3f,%.3f,%.3f], "depthpial":[%.3f,%.3f,%.3f], "dist":2.0, "soz":"0", "spikey":"0", "anat":"%s", "gridid":"%s", "PICS":"NaN"}]\n', subj, shortname, hem, elecmatrix(i,1), elecmatrix(i,2), elecmatrix(i,3), elecmatrix(i,1), elecmatrix(i,2), elecmatrix(i,3), anat, devname);
    else
    json = fprintf(fid,'{"subid":"%s", "elecid":"%s", "ElecType":"D", "Hem":"%s", "lepto":[%.3f,%.3f,%.3f], "depthpial":[%.3f,%.3f,%.3f], "dist":2.0, "soz":"0", "spikey":"0", "anat":"%s", "gridid":"%s", "PICS":"NaN"},', subj, shortname, hem, elecmatrix(i,1), elecmatrix(i,2), elecmatrix(i,3), elecmatrix(i,1), elecmatrix(i,2), elecmatrix(i,3), anat, devname);
    end
end
fclose(fid);