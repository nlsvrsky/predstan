function fitRandFunctions

ncores = str2num(getenv("NSLOTS"));
pool = parpool(ncores);

t = -3000:5:0;

fitParams = nan(4,100);

r1s = nan(100,length(t));
d1s = nan(100,length(t));
s1s = nan(100,length(t));
f1s = nan(100,length(t));

parfor ii=1:100
    tempOut = load(sprintf('output/randSequence/rand_out_%03d.mat',ii),'out');

    r1_rc = sum(tempOut.r1(:,end).*tempOut.stimList)./sum(tempOut.stimList)-mean(tempOut.r1(:,end));
    
    diffOfGamma = @(x) x(4) .* (t.*exp(t./x(1))-x(3).*t.*exp(t./x(2)));
    funcSSE = @(x) sum((diffOfGamma(x)-r1_rc).^2);

    funcFit = nan(4,100);
    funcVal = nan(1,100);
    
    for rr=1:100
        [tempFit,tempVal] = fminsearch(funcSSE,[randi(900),randi(900),1,1], ...
            optimset('MaxIter',1e4,'MaxFunEvals',1e4,'TolX',1e-12,'TolFun',1e-12));
        funcFit(:,rr) = tempFit;
        funcVal(rr) = tempVal;
    end

    bestFit = funcFit(:,find(funcVal==min(funcVal),1));
    fitParams(:,ii) = bestFit;

    r1s(ii,:) = r1_rc;
    d1s(ii,:) = sum(tempOut.d1(:,end).*tempOut.stimList)./sum(tempOut.stimList)-mean(tempOut.d1(:,end));
    s1s(ii,:) = sum(tempOut.s1(:,end).*tempOut.stimList)./sum(tempOut.stimList)-mean(tempOut.s1(:,end));
    f1s(ii,:) = sum(tempOut.f1(:,end).*tempOut.stimList)./sum(tempOut.stimList)-mean(tempOut.f1(:,end));
end

out.fitParams = fitParams;
out.r1 = r1s;
out.d1 = d1s;
out.s1 = s1s;
out.f1 = f1s;

save('output/randSeqFits.mat','out');
