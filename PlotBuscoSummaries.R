#' Assess completeness scores from a multisample busco run
#'
#' Plot scores from BUSCO summary files using compile_busco_summaries output. 
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta ortholog busco
#' @param dataframe is the compile_busco_summaries output
#' @param CSC "Complete_Single_Copy" = CSC,  # Default=Light Blue
#' @param CDC "Complete_Duplicated" = CDC,   # Default=Steel Blue
#' @param Frag "Fragmented" = Frag,            # Default=Medium Purple
#' @param Miss "Missing" = Miss                # Default=Peach Puff
#' @return Returns a dataframe of the summary scores, using each directory as an input ID
#' @export
#' @examples
#' compile_busco_summaries()

# Main function to process all directories
PlotBuscoSummaries <- function(dataframe, CSC="#89CFF0", CDC="#4682B4", Frag="#9370DB", Miss="#FFCC99") {
# Prepare data for plotting
plot_data <- all_summaries %>%
  arrange(desc(Complete)) %>%  # Sort by Complete percentage
  pivot_longer(cols = c(Complete_Single_Copy, Complete_Duplicated, Fragmented, Missing),
               names_to = "Category", values_to = "Percentage") %>%
  mutate(
    Category = factor(Category, levels = c("Missing", "Fragmented", "Complete_Duplicated", "Complete_Single_Copy")),
    ID = factor(ID, levels = rev(unique(all_summaries$ID[order(-all_summaries$Complete)])))  # Reverse factor levels
  )

# Define custom colors
colors <- c(
  "Complete_Single_Copy" = CSC,  # Default=Light Blue
  "Complete_Duplicated" = CDC,   # Default=Steel Blue
  "Fragmented" = Frag,            # Default=Medium Purple
  "Missing" = Miss                # Default=Peach Puff
)

# Plot
ggplot(plot_data, aes(y = ID, x = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Percentage",
    y = "Sample ID",
    fill = "Category",
    title = "BUSCO Results by Sample"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
}

