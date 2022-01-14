excluded$reasons


reasons <- recode_tags(list(no_measure=no_measure, 
                            invalid_comparator=invalid_comparator, 
                            included=included,
                            incorrect_habitat=incorrect_habitat,
                            incorrect_geog=incorrect_geog,
                            invalid_time=invalid_time), excluded$reasons)

excluded$reasons[reasons==""]

no_measure <- c("no measure")
invalid_comparator <- c("comparing")
included <- c("included")
incorrect_habitat <- c("forest")
incorrect_geog <- c("North")
invalid_time <- c("winter")



library(DiagrammeR)


DiagrammeR::grViz("digraph {

graph [layout = dot, rankdir = LR]

# define the global styles of the nodes. We can override these in box if we wish
node [shape = rectangle, style = filled, fillcolor = Linen]

data1 [label = 'Dataset 1', shape = folder, fillcolor = Beige]
data2 [label = 'Dataset 2', shape = folder, fillcolor = Beige]
process [label =  'Process \n Data']
statistical [label = 'Statistical \n Analysis']
results [label= 'Results']

# edge definitions with the node IDs
{data1 data2}  -> process -> statistical -> results
}")