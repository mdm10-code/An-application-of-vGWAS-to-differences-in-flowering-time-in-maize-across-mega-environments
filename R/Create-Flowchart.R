##################################################################
#################### Create Flowchart ############################
##################################################################

#As always, load your libraries
library(DiagrammeR)

#This is an example R Visualizations: Flow Charts in R by Paul Aleksis

grViz(diagram = "digraph flowchart {
  node [fontname = arial, shape = oval]
  tab1 [label = '@@1']
  tab2 [label = '@@2']
  tab3 [label = '@@3']
  tab4 [label = '@@4']
  
  tab1 -> tab2 -> tab3 -> tab4;
}
  
  [1]: 'Learning Data Science'
  [2]: 'Industry vs Technical Knowledge'    
  [3]: 'Statistics vs Mathematics Knowledge' 
  [4]: 'Brown-Forsythe Test'
  ")

##This example is for branching
pdf("Flowchart-for-Analyses.pdf", width = 10)
grViz(diagram = "digraph flowchart {
      # define node aesthetics
      node [fontname = Arial, shape = oval, color = Lavender, style = filled]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
# set up node layout
      tab1 -> tab2;
      tab2 -> tab3;
      tab2 -> tab4
      }
[1]: 'Learning Data Science'
      [2]: 'Industry vs Technical Knowledge'
      [3]: 'Python/R'
      [4]: 'Domain Experience'
      ")

###With color


pdf("Flowchart-for-Analyses.pdf", width = 10)
grViz(diagram = "digraph flowchart {
        # graph attributes
  graph [overlap = true]
      
      # define node aesthetics
      node [fontname = Helvetica, shape = rectangle, color = darkorange3, style = filled, fontcolor = blue4]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      tab6 [label = '@@6']
      tab7 [label = '@@7']
      tab8 [label = '@@8']
# set up node layout
      tab1 -> tab2;
      tab2 -> tab3;
      tab3 -> tab4
      tab4 -> tab5
      tab5 -> tab6
      tab6 -> tab7
      tab6 -> tab8
      }
      [1]: 'Calculate Best Linear Unbiased Predictors (BLUPs) for each mega-environment'
      [2]: 'Take the difference of the BLUPs between mega-environments for the genotype-by-environment (GxE) response variable'
      [3]: 'Run Principal Component Analysis (PCAs) and Van Raden Kinship Matrix in TASSEL'
      [4]: 'Run BIC Criterion to determine number of PCAs to use for population structure'
      [5]: 'Run the Unified Mixed Linear Model in TASSEL'
      [6]: 'Extract residuals from the Unified Mixed Linear Model'
      [7]: 'Run the Brown-Forsythe Test (BFT)'
      [8]: 'Run the Double Generalized Linear Model (DGLM)'
      ")
dev.off()


###Revision #2
pdf("Flowchart-for-Analyses.pdf", width = 6)
grViz(diagram = "digraph flowchart {
        # graph attributes
  graph [overlap = true]
      
      # define node aesthetics
      node [fontname = Helvetica, shape = rectangle, color = darkorange3, style = filled, fontcolor = blue4]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      tab6 [label = '@@6']
      tab7 [label = '@@7']
      tab8 [label = '@@8']
# set up node layout
      tab1 -> tab2;
      tab2 -> tab4;
      tab3 -> tab4
      tab4 -> tab5
      tab5 -> tab6
      tab6 -> tab7
      tab6 -> tab8
      }
      [1]: 'Calculate Best Linear Unbiased Predictors (BLUPs) for each mega-environment'
      [2]: 'Take the difference of the BLUPs between mega-environments for the genotype-by-environment (GxE) response variable'
      [3]: 'Run Principal Component Analysis (PCA) and VanRaden Kinship Matrix in TASSEL'
      [4]: 'Use BIC to determine number of principal components to account for population structure'
      [5]: 'Run the Unified Mixed Linear Model in TASSEL'
      [6]: 'Extract residuals from the Unified Mixed Linear Model'
      [7]: 'Run the Brown-Forsythe Test (BFT)'
      [8]: 'Run the Double Generalized Linear Model (DGLM)'
      ")
dev.off()

