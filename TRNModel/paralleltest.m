cluster = parcluster
job = createCommunicatingJob(cluster,'Type','pool');
% inputs = {trnmodel,FBAmodel,Amat,Bmat,Cmat};
createTask(job,@rand,1,{3});
submit(job);
% wait(job);
% data = fetchOutputs(job);