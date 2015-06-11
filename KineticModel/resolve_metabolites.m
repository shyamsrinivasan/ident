%resolve_metabolites(FBAmodel)
load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Kinetic Model\centralIshiiModel.mat');
load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Kinetic Model\toy_run_1.mat');
mets = lower(model.mets);
for ime = 1:length(mets)
    mets{ime} = [mets{ime} '[c]'];
end
MC = zeros(length(FBAmodel.Metabolites),1000);
for im = 1:length(FBAmodel.Metabolites)
    tfm = strcmpi(FBAmodel.Metabolites{im},mets);
    if any(tfm)
        fprintf('%d Metabolite:%s\n',im,FBAmodel.Metabolites{im});
        MC(im,:) = points(tfm,:);
    end
    iv = iv + im;
end