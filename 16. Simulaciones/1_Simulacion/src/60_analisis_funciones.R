# a %>% 
#   select(group, PTA) %>% 
#   unnest(PTA) %>% 
#   ggplot(aes(x = MIC, y = PTA, col = group)) +
#   geom_line() + 
#   scale_x_continuous(
#     trans = 'log2',
#     breaks = trans_breaks("log2", function(x) 2^x),
#     labels = trans_format("log2", math_format(2^.x)))
# 
# a %>% 
#   select(group, fTmasMIC) %>% 
#   unnest(fTmasMIC) %>% 
#   pivot_longer(cols = matches('^ID'),names_to = 'ID', values_to = 'ind_pkpd') %>% 
#   ggplot(aes(x = MIC, y = ind_pkpd, col = ID)) +
#   geom_line() + 
#   scale_x_continuous(
#     trans = 'log2',
#     breaks = trans_breaks("log2", function(x) 2^x),
#     labels = trans_format("log2", math_format(2^.x))) +
#   facet_wrap(.~group, ncol = 2) +
#   theme(legend.position = "none")
