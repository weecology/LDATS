library(tidyverse)
library(ggplot2)

N = 10^2
N_topics = 3
N_docs = 2
N_species = 3
alpha = 2
alpha_beta = .2


set.seed(21)
beta = matrix(c(0.64, 0.36, 0, 
                0, 0.66, 0.34, 
                0, .14, .86),
              N_topics, N_species, byrow = TRUE)

set.seed(1)
thetas = list(c(.11, .45, .44),
              c(.24, .46, .30))

zs = thetas %>% 
  map(~sample.int(n = N_topics, size = N, prob = .x, replace = TRUE))


pick_species = function(z){
  map_dbl(z, ~sample.int(n = N_species, size = 1, prob = beta[.x, ]))
}

assemblages = map(zs, pick_species)

get_xmin = function(doc){
  (doc - 1) * (sqrt(N) + 1)
}

labels = round(beta * 100) %>% 
  apply(1, list) %>% 
  flatten() %>% 
  map(~paste(.x, "%  ", sep = "")) %>% 
  map(paste, collapse = "\n") %>% 
  flatten_chr()

output = transpose(.l = list(doc = 1:N_docs,
                             z = zs, 
                             word = assemblages)) %>% 
  map(as_data_frame) %>% 
  bind_rows() %>% 
  arrange(doc, z, word) %>% 
  group_by(doc) %>% 
  mutate(x = get_xmin(doc) + rep(1:sqrt(N), sqrt(N)), 
         y = rep(sqrt(N):1, each = sqrt(N)))

rects = data_frame(
  xmin = rep(get_xmin(seq_len(N_docs)) + .5, 2),
  xmax = rep(get_xmin(seq_len(N_docs)) + sqrt(N) + .5, 2),
  ymin = rep(c(.5, 1 + sqrt(N)), each = 2),
  ymax = rep(c(.5 + sqrt(N), 3 + sqrt(N)), each = 2)
)

topic_positions = seq(min(rects), max(rects), length = N_topics + 2)
topic_positions = topic_positions[-c(1, length(topic_positions))]

arrows = data_frame(
  x = rep(topic_positions, times = N_docs),
  y = sqrt(N) + 9,
  xend = rep(unique(rects$xmin + rects$xmax), each = N_topics) / 2 + 
    1:N_topics - mean(1:N_topics),
  yend = .25 + max(rects$ymax),
  size = 5 * flatten_dbl(thetas)
)

counts = output %>% 
  group_by(doc, word) %>% 
  summarize(n = n()) %>% 
  group_by(doc) %>% 
  summarize(counts = paste0(n, collapse = "\n"))

left_labels = tribble(
  ~label,                  ~y,
  "Topic\ndefinitions",    max(arrows$y) + 2,
  "Topic\nproportions",    min(arrows$yend),
  "Assemblages",           sqrt(N)/2 + 1,
  "Species\ncomposition",   -1.5
)
  

ggplot() +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, 
                              ymin = ymin, ymax = ymax),
            color = "black", fill = NA) + 
  geom_point(data = output, aes(x = x,
                                y = y, 
                                color = factor(z), 
                                shape = factor(word)), 
             size = 6) + 
  scale_shape_discrete("Species", labels = LETTERS[1:N_species]) + 
  scale_color_brewer(palette = "Dark2", guide = FALSE) +
  scale_fill_brewer("Community type", palette = "Dark2") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.background = element_blank(), panel.spacing = unit(2, "line")) +
  coord_equal(xlim = c(-6, (sqrt(N) + 1) * N_docs),
              ylim = c(-4, max(arrows$y) + sqrt(N) / 2)) + 
  xlab("") +
  ylab("") +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend,
                                  size = size), 
               arrow = arrow(type = "closed"),
               size = arrows$size) +
  geom_tile(data = arrows[1:N_topics, ], 
            aes(x = x, y = y + 1, height = sqrt(N) / 2, width = sqrt(N) / 2, 
                fill = factor(1:N_topics))) +
  geom_text(data = arrows[1:N_topics, ], 
             aes(x = x, y = y + 1, label = labels),
             size = 8) + 
  geom_text(data = left_labels, aes(x = 0, y = y, label = label), 
            size = 8, hjust = "right", vjust = "top") + 
  geom_text(aes(x = sqrt(N) * seq(0.5, N_docs - .5) + seq_len(N_docs), 
                y = -1, label = counts$counts), 
            size = 8, vjust = 1, hjust = 0)

