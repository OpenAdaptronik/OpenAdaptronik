m = 1;
k = 5;
d = .5;


args = sciargs();
iarg = find(args=='-args');
m = strtod(args(iarg+1));
k = strtod(args(iarg+2));
d = strtod(args(iarg+3));


exec('transff.sce');
[w0 wd]=transff(m, k, d);
