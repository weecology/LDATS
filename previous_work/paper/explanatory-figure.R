library(tidyverse)
library(ggplot2)

N = 10^2
N_topics = 3
N_docs = 2
N_species = 3
alpha = 2
alpha_beta = .2

# From the 4-class PuOr palette on Colorbrewer2.org
colors = c('#e66101','#fdb863','#b2abd2','#5e3c99')[c(3, 1, 2)]
colors = viridis::viridis(12)[c(7, 11, 9)]


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
  xmin = get_xmin(seq_len(N_docs)) + .5,
  xmax = get_xmin(seq_len(N_docs)) + sqrt(N) + .5,
  ymin = rep(.5, 2),
  ymax = rep(.5 + sqrt(N), 2),
  fill = factor(1, levels = 1:N_topics),
  alpha = 0)

proportion_rects = lapply(1:N_docs,
                          function(i) {
                            data_frame(xmin = rects$xmin[i] + 
                                         c(0, thetas[[i]][-N_topics]) * sqrt(N),
                                       xmax = c(xmin[-1], rects$xmax[[i]]),
                                       ymin = sqrt(N) + 1.5,
                                       ymax = sqrt(N) + 2.5,
                                       fill = factor(1:N_topics),
                                       alpha = 1)
                          }) %>% 
  bind_rows()


topic_positions = seq(min(rects$xmin), max(rects$xmax), length = N_topics + 2)
topic_positions = topic_positions[-c(1, length(topic_positions))]


arrows = data_frame(
  x = rep(topic_positions, times = N_docs),
  y = max(proportion_rects$ymax) + 7,
  xend = (proportion_rects$xmin + proportion_rects$xmax) / 2,
  yend = .1 + max(proportion_rects$ymax),
  size = 1
)

counts = output %>% 
  group_by(doc, word) %>% 
  summarize(n = n()) %>% 
  group_by(doc) %>% 
  summarize(counts = paste0(n, collapse = "\n"))

composition_x = sqrt(N) * seq(0.5, N_docs - .5) + seq_len(N_docs)

left_labels = tribble(
  ~label,                  ~y,
  "A  Community-\ntype definitions",    max(arrows$y) + 2,
  "B  Community-\ntype proportions",    min(arrows$yend) + 0.5,
  "C  Assemblages",           sqrt(N)/2 + 1,
  "D  Species\ncomposition",   -1.5
)
point_size = 3.4

basic_plot = ggplot() +
  geom_segment(data = arrows, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               arrow = arrow(type = "closed", length = unit(.1, "inches")),
               size = 1) + 
  geom_rect(data = bind_rows(rects, proportion_rects), 
            aes(xmin = xmin, xmax = xmax, 
                ymin = ymin, ymax = ymax, fill = fill, alpha = alpha),
            color = "black") + 
  geom_tile(data = arrows[1:N_topics, ], 
            aes(x = x, y = y + 1, height = sqrt(N) / 2, width = sqrt(N) / 2, 
                fill = factor(1:N_topics))) +
  geom_point(data = output, 
             aes(x = x,
                 y = y, 
                 color = factor(z), 
                 shape = factor(word)), 
             size = point_size) + 
  geom_point(aes(x = rep(composition_x, each = N_species),
                 shape = rep(factor(1:N_species), N_docs), 
                 y = rep(0.25 - 1.6 * (1:N_species), N_docs)),
             size = point_size) + 
  geom_point(aes(x = rep(unique(arrows$x) + 1.4, each = N_topics), 
                 y = rep(max(arrows$y) + 4.2 - 1.6 * (1:N_species), N_topics),
                 shape = rep(factor(1:N_species), N_topics)),
             size = point_size)

styled_plot = basic_plot + 
  scale_shape_discrete("Species", labels = paste("Sp.", 1:N_species)) + 
  scale_color_manual(values = colors, guide = FALSE) +
  scale_alpha_continuous(range = c(0, 1), guide = FALSE) +
  scale_fill_manual("Community type", values = colors) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.background = element_blank()) +
  coord_equal(xlim = c(-7, (sqrt(N) + 1) * N_docs),
              ylim = c(-4, max(arrows$y) + sqrt(N) / 2)) + 
  xlab("") +
  ylab("") +
  theme(legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12))

final_plot = styled_plot + 
  geom_text(data = arrows[1:N_topics, ], 
            aes(x = x, y = y + 1, label = labels),
            size = 4.5) + 
  geom_text(data = left_labels, aes(x = 0, y = y, label = label), 
            size = 4.5, hjust = "right", vjust = "top") + 
  geom_text(aes(x = composition_x - 1, 
                y = -1, label = counts$counts), 
            size = 4.5, vjust = 1, hjust = 1)


final_plot
ggsave("explanatory.tiff", width = 8.5, height = 6)
