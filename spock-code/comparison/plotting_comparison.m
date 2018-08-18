% High level strategy
% 1. Loop through all mat files (each experiment replicate)
% 2. Collect data points and gene names
% 3. Fill in missing data points from replicates 1 and 2 with anything
%    present in replicate 3
% 4. Save to excel file, save additional sheets with standard dev. removed

close all; clear all;
folders = {'Matfiles2','Matfiles3','Matfiles4'};
filename2 = 'gene_final.xlsx';


plotting = [];
gene = {};
expM = zeros(1,2);
expS = expM;
for z = 1:length(folders) %%loop through folders
    load(['colavg',num2str(z)]);
    genecolIndiv = {};
    
    fileList = dir(folders{z});
    data_pt = [];
    for k = 1:length(fileList)
        
        % check if current item is a directory
        if fileList(k).isdir
            continue
        elseif strcmp(fileList(k).name,'.DS_Store')
            continue
        end
        % check if it's a PC or Mac to build file structure
        if ispc
            str = load([fileList(k).folder,'\',fileList(k).name]);
        else
            str = load([fileList(k).folder,'/',fileList(k).name]);
        end
        
        %dynamically load struct from file
        structName = fieldnames(str);
        structName = structName{1};
        experiment = str.(structName);
        
        
        grid = experiment.grid;
        
        
        disp('Colony Filtering')
        for r = 1:size(grid,1)
            for c = 1:size(grid,2)
                if isempty(experiment.grid(r,c).mean_colony_intensity)
                    experiment.grid(r,c).mean_colony_intensity = 0;
                end
                data_pt(end+1) = experiment.grid(r,c).mean_colony_intensity./mean(colavg);
                genecolIndiv{end+1} = [experiment.grid(r,c).geneName,experiment.grid(r,c).geneLoc,experiment.grid(r,c).genePlate];

            end
        end
        
    end
    expM(z) = mean(data_pt);
    expS(z) = std(data_pt);
    plotting = [plotting;data_pt];
    gene{z} = genecolIndiv;
    
end


% define gene name, plate, and well location to write to output table
geneall = gene{1};
gene_names = cell(1,length(geneall));
gene_well = cell(1,length(geneall));
gene_plate = cell(1,length(geneall));
for k = 1:length(geneall)
    gene_names{k} = geneall{k}{1};
    gene_well{k} = geneall{k}{2};
    gene_plate{k} = geneall{k}{3};
end
output_table1 = table(plotting(1,:)',plotting(2,:)',plotting(3,:)',gene_names',gene_well',gene_plate');
output_table1.Properties.VariableNames = {'rep1','rep2','rep3','gene_name','gene_well','gene_plate'};
writetable(output_table1,filename2,'Sheet',1)

% find any index from replicate # 3 that is not equal to 0, where there
% exists a 0 in either replicates 1 or 2
ind1 = find(plotting(1,:)==0);
ind2 = find(plotting(2,:)==0);
ind3 = find(plotting(3,:)~=0);


% fill in missing points for experiment 2
inter1 = intersect(ind3,ind2);
p3 = zeros(3,size(plotting,2));
p3(3,inter1) = plotting(3,inter1);
plotting(2,inter1) = plotting(2,inter1)+plotting(3,inter1);
% fill in missing points for experiment 1
inter2 = intersect(ind3,ind1);
p3 = zeros(3,size(plotting,2));
p3(3,inter2) = plotting(3,inter2);
plotting(1,inter2) = plotting(1,inter2)+plotting(3,inter2);

new_plotting = zeros(2,length(plotting));
for k = 1:length(plotting)
    
    diff1 = abs(plotting(1,k)-plotting(2,k));
    diff2 = abs(plotting(2,k)-plotting(3,k));
    diff3 = abs(plotting(1,k)-plotting(3,k));
    
    [~,i] = min([diff1, diff2, diff3]);
    
    switch i
        case 1
            new_plotting(:,k) = plotting(1:2,k);
        case 2
            new_plotting(:,k) = plotting(2:3,k);
        case 3
            new_plotting(:,k) = plotting([1,3],k);
    end
    
end

plotting = new_plotting;
output_table2 = table(plotting(1,:)',plotting(2,:)',gene_names',gene_well',gene_plate');
output_table2.Properties.VariableNames = {'rep1','rep2','gene_name','gene_well','gene_plate'};
writetable(output_table2,filename2,'Sheet',2)

figure; hold on, scatter(plotting(1,:),plotting(2,:)),set(gca,'fontsize',30),xlabel('Rep 1'),ylabel('Rep 2')
xlim([0 3]), ylim([0 3]), axis square
expM = [mean(plotting(1,:)),mean(plotting(1,:))]; expS = [std(plotting(1,:)),std(plotting(2,:))];


% remove any remaining zeros
ind1 = find(plotting(1,:)==0);
ind2 = find(plotting(2,:)==0);
ind = unique([ind1,ind2]);
plotting(:,ind) = [];
geneall(ind) = [];



% initiate cell
meani = mean(plotting);
avg = {};
for q = 1:length(meani)
    avg{end+1} = meani(q);
end

% removing points greater than 2 standard deviations, save to separate sheet
ind1 = find(plotting(1,:)> (expM(1) + 2.*expS(1)));
ind2 = find(plotting(2,:)> (expM(2) + 2.*expS(2)));
ind = union(ind1,ind2);
ind = [ind, setdiff(ind,ind3)];
filename = 'colonydata.xlsx';
if ~isempty(geneall(ind))
    T1 = cell2table([avg(ind)',geneall(ind)']);
    writetable(T1,filename2,'Sheet',3)
end

% removing points less than 2 standard deviations, save to separate sheet
ind1 = find(plotting(1,:)< (expM(1) - 2.*expS(1)));
ind2 = find(plotting(2,:)< (expM(2) - 2.*expS(2)));
ind = [];
for q = 1:length(ind2)
    for qq = 1:length(ind1)
        if ind1(qq) == ind2(q)
            ind(end+1) = ind1(qq);
        end
    end
end
if ~isempty(geneall(ind))
    T2 = cell2table([avg(ind)',geneall(ind)']);
    writetable(T2,filename2,'Sheet',4)
end


figure; plot(plotting(1,:),'o'),xlim([0 4100]),ylim([0 4]), set(gca,'fontsize',30),title('Rep 1'),xlabel('Index'),ylabel('Avg intensity')%,plotting(2,:))
figure; plot(plotting(2,:),'o'),xlim([0 4100]),ylim([0 4]), set(gca,'fontsize',30),title('Rep 2'),xlabel('Index'),ylabel('Avg intensity')%,plotting(2,:))

