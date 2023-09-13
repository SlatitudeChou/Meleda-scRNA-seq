### cell portion test
devtools::install_github("rpolicastro/scProportionTest")
library("scProportionTest")

seurat_data <- system.file("extdata", "example_data.RDS", package = "scProportionTest")
seurat_data <- readRDS(seurat_data)

seurat_data@meta.data %>% head()
seurat_data <- seurat_immune

prop_test <- sc_utils(seurat_data)
prop_test_AD <- permutation_test(
  prop_test, cluster_identity = "typesubfine",
  sample_1 = "ADL", sample_2 = "ADHC",
  sample_identity = "group"
)
p1 <- permutation_plot(prop_test_AD) + NoLegend()

prop_test_AS <- permutation_test(
  prop_test, cluster_identity = "typesubfine",
  sample_1 = "ASL", sample_2 = "ASHC",
  sample_identity = "group"
)
p2 <- permutation_plot(prop_test_AS) + NoLegend()

prop_test_CRS <- permutation_test(
  prop_test, cluster_identity = "typesubfine",
  sample_1 = "CRSL", sample_2 = "CRSHC",
  sample_identity = "group"
)
p3 <- permutation_plot(prop_test_CRS) + NoLegend()

plot_grid(p1,p2,p3,nrow = 1)

ggsave("./results_immune/portiontest.jpeg", device = 'jpeg', width = 7, height = 1.8, unit = "in", dpi = 300, limitsize = FALSE)


library(ggplot2)

# create a dataset
group <- c(rep("ADHC" , 7) ,rep("ADL" , 7),
           rep("ASHC" , 7), rep("ASL" , 7),
           rep("CRSHC" , 7), rep("CRSL" , 7))
main <- rep(c("B", 'CD4T','CD8T/NKT','DC','Macro' ,'Mast', 'Proliferation') , 6)
number <- as.vector(table(seurat_immune$group_subfine))
number <- c(number[1:34],c(0),number[35:41])
number <- number[c(-7,-14,-21,-28,-41)]
immuneportion <- data.frame(group, main, number)
sum(number)


# Stacked + percent
ggplot(immuneportion, aes(fill=main, y=number, x=group)) + 
  geom_bar(position="fill", stat="identity") + 
  theme_set(theme_bw()) 
ggsave("./results_immune/portionplot.jpeg", device = 'jpeg', width = 5, height = 5, unit = "in", dpi = 300, limitsize = FALSE)

immunenumber <- data.frame(group = names(table(seurat_immune$group)), 
                           number = as.vector(table(seurat_immune$group)))

immunenumber$numberlog <- log(immunenumber$number)
  
ggplot(immunenumber, aes(x=group, y=number, fill=group)) +
  geom_bar(stat="identity")+theme_minimal()+ 
  theme_set(theme_bw()) 
ggsave("./results_immune/numberplot.jpeg", device = 'jpeg', width = 5, height = 5, unit = "in", dpi = 300, limitsize = FALSE)


