#pseudogenes_file = '\\\\wsl$\\Ubuntu\\home\\hoeps\\turbo_mount_3\\data\\research\\projects\\wolfram\\projects\\paraphase_find_indels\\final_table.txt'

downsample_results = '\\\\wsl$\\Ubuntu\\home\\hoeps\\turbo_mount_3\\data\\research\\projects\\wolfram\\projects\\paraphase_find_indels\\final_table.txt'

data = read.table(pseudogenes_file, header=1, sep=' ')

library(ggplot2)
library(reshape2)
library(reshape)

dm = melt(data)


dm[dm$value>10,'value'] = 10
ggplot(dm) + geom_boxplot(aes(x=Metric, y=value, fill=Metric)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_continuous(breaks=seq(0,10,2))

# Sort boxplot by median
ggplot(dm) + geom_boxplot(aes(x=reorder(Metric, value, median), y=value, fill=Metric)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_continuous(breaks=seq(0,10,2))