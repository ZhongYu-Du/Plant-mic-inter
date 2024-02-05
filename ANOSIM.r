library(vegan)


dis <- read.delim('bray.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
dis <- as.dist(dis)	
otu <- read.delim('otu_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))

group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE)


anosim_result_dis <- anosim(dis, group$site, permutations = 999)	sult_otu <- anosim(otu, group$site, permutations = 999, distance = 'bray') 	#?????????eosim_result_dis)
#????
nameim_result_dis)
anosim_result_dis$signif	#p ֵ
anosim_result_dis$statistic	#R ֵ

#??ͼչʾ
(group$site)

dir.create('anosim_two', recursive = TRUE)
anosim_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
	for (j in (i + 1):length(group_name)) {
		group_ij <- subset(group, site %in% c(group_name[i], group_name[j]))
		otu_ij <- otu[group_ij$names, ]
		anosim_result_otu_ij <- anosim(otu_ij, group_ij$site, permutations = 999, distance = 'bray')	#?????û????? 999 ??R ֵ??u_ij$signif <= 0.001) Sig <- '***'
		else if (anosim_result_otu_ij$signif <= 0.01) Sig <- '**'
		else if (anosim_result_otu_ij$signif <= 0.05) Sig <- '*'
		else Sig <- NA
		anosim_result_two <- rbind(anosim_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 'Bray-Curtis', anosim_result_otu_ij$statistic, anosim_result_otu_ij$signif, Sig))
		
		#ÿ??ѭ??????ͼ? p ֵ? data.frame(anosim_result_two, stringsAsFactors = FALSE)
names(anosim_result_two) <- c('group', 'distance', 'R', 'P_value', 'Sig')
write.table(anosim_result_two, 'anosim_two/ANOSIM.ret', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
