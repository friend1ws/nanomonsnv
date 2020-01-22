library(tidyverse)

D <- read_tsv("~/Downloads/out3.mini.txt", col_names = FALSE, guess_max = Inf) %>% filter(X11 >= -log10(0.01), X19 != "Ambiguous")


D$X15[D$X15 == "---"] <- "0"
D$X16[D$X16 == "---"] <- "0"

D$X15_X16 <- pmin(as.numeric(D$X15), as.numeric(D$X16))

D1 <- D %>% group_by(X15_X16, X19) %>% summarize(Count = n())


ggplot(D, aes(x = X11, fill = X19)) + 
  geom_histogram(position = "dodge", binwidth = 0.25) +
  theme_bw() +
  labs(fill = "", x = "log10 (P-value of the Fisher's exact test)") +
  theme(legend.position = "bottom") +
  xlim(c(2, 12))


ggplot(D1 %>% filter(X15_X16 >= 1), aes(x = X15_X16, y = Count, fill = X19)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(fill = "", x = "Min (Max (Positive strand base qualities), Max(Negative strand base qualities))") +
  theme(legend.position = "bottom") +
  xlim(c(1, 32))

ggplot(D, aes(x = X12, fill = X19)) + 
  geom_histogram(position = "dodge", binwidth = 0.025) +
  theme_bw() +
  labs(fill = "", x = "The ratio of positive strand variant reads") +
  theme(legend.position = "bottom") 

ggplot(D %>%, aes(x = X17, fill = X19)) + 
  geom_histogram(position = "dodge", binwidth = 0.025) +
  theme_bw() +
  labs(fill = "", x = "Maximum of the non-matched control error ratios") +
  theme(legend.position = "bottom") +
  xlim(c(0, 0.5))


ggplot(D %>%, aes(x = X18, fill = X19)) + 
  geom_histogram(position = "dodge", binwidth = 0.025) +
  theme_bw() +
  labs(fill = "", x = "Median of the non-matched control error ratios") +
  theme(legend.position = "bottom") +
  xlim(c(0, 0.5))


D_filt <- D %>% filter(X12 <= 0.95, X12 >= 0.05,X15_X16 >= 10, X17 <= 0.12, X18 <= 0.1)


ggplot(D_filt, aes(x = X11, fill = X19)) + 
  geom_histogram(position = "dodge", binwidth = 0.25) +
  theme_bw() +
  labs(fill = "", x = "log10 (P-value of the Fisher's exact test)") +
  theme(legend.position = "bottom") +
  xlim(c(2, 16))



