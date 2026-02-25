# Make sure you have the package installed
# install.packages("DiagrammeR")

library(DiagrammeR)

# Create the DAG using the grViz function and DOT language
goby_dag_viz <- grViz("
digraph goby_sem {
  # --- Graph Attributes ---
  graph [layout = dot, rankdir = LR, splines=true]
  
  # --- Node Definitions & Styling ---
  # Default node style
  node [shape = box, style = rounded, fontname = Helvetica]

  # Subgraph for Latent Variables
  subgraph cluster_latent {
    label = 'Latent Drivers';
    style = invis;
    U_Physical [label = 'U (Physical)', fillcolor = grey88, style=filled];
    U_Biological [label = 'U (Biological)', fillcolor = grey88, style=filled];
  }

  # Subgraph for Exogenous Predictors
  subgraph cluster_exogenous {
    label = '';
    style = invis;
    # Environmental drivers
    node [fillcolor = '#5bc0de', style=filled, fontcolor=white];
    Rain; Wind; Substrate; Goby_lag;
    
    # NEW: Temporal & Micro-habitat drivers
    node [fillcolor = '#9b59b6', style=filled, fontcolor=white];
    Micro; Year;
  }

  # Subgraph for Mediators
  subgraph cluster_mediators {
    label = '';
    style = invis;
    # Physical Mediators
    node [fillcolor = '#f0ad4e', style=filled, fontcolor=white];
    BreachDays; Temp; DO;
    
    # Biological Mediators
    node [fillcolor = '#5cb85c', style=filled, fontcolor=white];
    SAV; SC_count; SB_count;
  }

  # Outcome Node
  Goby [fillcolor = '#d9534f', style=filled, fontcolor=white];

  # --- Edge (Arrow) Definitions ---
  # Correlated Error
  Rain -> Wind [dir=both, color=grey50, style=dashed];

  # Latent Variable Paths
  U_Physical -> {BreachDays, Temp, DO}
  U_Biological -> {SAV, SC_count, SB_count}
  
  # Causal Cascade
  Rain -> BreachDays -> Temp -> DO -> SAV;
  Wind -> {Temp, DO};
  Temp -> SAV;
  
  # Biological Sub-models
  {DO, SAV, Substrate} -> SC_count;
  {DO, SAV} -> SB_count;
  
  # Paths to the final Goby outcome
  {BreachDays, Temp, DO, SAV, SC_count, SB_count, Substrate, Micro, Year, Goby_lag} -> Goby;
}
")

# Print the DAG
goby_dag_viz




########--- updated 2025-12-28
# DAG with Wind/Rain at top, cascading down to Goby at bottom
library(dagitty)
library(ggdag)
library(ggplot2)

dag_string <- '
dag {
  Goby [outcome,pos="3,0"]
  SAV [pos="2,2"]
  DO [pos="3,3"]
  Temp [pos="2,4"]
  BreachDays [pos="1,5"]
  SC [pos="1.5,1.5"]
  SB [pos="4,1.5"]
  Wind [pos="4,6"]
  Rain [pos="0,7"]
  U_bio [latent,pos="5,2"]
  U_phys [latent,pos="5,5"]
  
  Rain -> BreachDays
  BreachDays -> Temp
  Wind -> Temp
  Wind -> DO
  Temp -> DO
  DO -> SAV
  Temp -> SAV
  SAV -> SC
  DO -> SC
  SAV -> SB
  DO -> SB
  
  SAV -> Goby
  DO -> Goby
  Temp -> Goby
  BreachDays -> Goby
  SC -> Goby
  SB -> Goby
  
  U_bio -> SAV
  U_bio -> SC
  U_bio -> SB
  U_phys -> Temp
  U_phys -> DO
  U_phys -> BreachDays
  U_phys -> Wind
  U_phys -> Rain
}
'

dag_dagitty <- dagitty(dag_string)

# Base plot
plot(dag_dagitty, main = "SEM DAG: Wind/Rain → Goby")
