import pandas as pd
import glob
import os
import numpy as np
from sklearn.preprocessing import QuantileTransformer
import matplotlib.pyplot as plt

mut_list = pd.read_excel('selected_mutations.xlsx')

file_list = glob.glob('/your-dir/*.txt') ### List of IR result files "IRFinder-IR-nondir.txt" named after the samples.

# Function for reading multiple files and store them seperately based on intron depth information and IR ratios.
def data_reader(fil_list,ir_locs,depth_locs):
  ir_df = pd.DataFrame()
  depth_df = pd.DataFrame()
  # use a loop to read each file and append specific columns to the final dataframe
  for i, file in enumerate(file_list):
      # read in the file as a dataframe
      df = pd.read_csv(file, sep='\t')
      #if 'SRR' not in file:
        # if it is the first file, select the first four columns and the twentieth column
      if i == 0:
              ir_df = pd.concat([ir_df, df.iloc[:, ir_locs]], axis=1)
              depth_df = pd.concat([depth_df, df.iloc[:,depth_locs]], axis=1)
              # get the name of the file before .csv and use it to rename the twentieth column
              column_name = os.path.splitext(os.path.basename(file))[0]
              ir_df = ir_df.rename(columns={ir_df.columns[-1]: column_name})
              depth_df = depth_df.rename(columns={depth_df.columns[-1]: column_name})
          # for subsequent files, select only the twentieth column
      else:
              new_column_ir = pd.DataFrame(df.iloc[:, 19])
              new_column_depth = pd.DataFrame(df.iloc[:, 8])
              # get the name of the file before .csv and use it to rename the new column
              column_name = os.path.splitext(os.path.basename(file))[0]
              new_column_ir = new_column_ir.rename(columns={new_column_ir.columns[-1]: column_name})
              new_column_depth = new_column_depth.rename(columns={new_column_depth.columns[-1]: column_name})
              ir_df = pd.concat([ir_df, new_column_ir], axis=1)
              depth_df = pd.concat([depth_df, new_column_depth], axis=1)

  # display the final dataframe
  return ir_df, depth_df
  
  
# Function to calculate CPM for a single sample
def calculate_cpm(sample_counts):
    library_size = sample_counts.sum()
    cpm = (sample_counts / library_size) * 1_000_000
    return cpm
    
    
# Column numbers of the IR (19) and depth (8) information combination of the meta-data columns (0,1,2,3)    
ir_locs = [0, 1, 2, 3, 19]
depth_locs =  [0, 1, 2, 3, 8]

ir_df, depth_df = data_reader(file_list,ir_locs,depth_locs)

# Separate the first four non-integer columns from the remaining ones
metadata_columns = ir_df.iloc[:, :4]
ir_data = ir_df.iloc[:, 4:]

# Assign weights to each column (sample) based on the column names
# Define the default_weight and higher_weight
default_weight = 1.0
higher_weight = 2.0

weights = np.array([higher_weight if col.startswith('SRR') else default_weight for col in ir_data.columns])

# Calculate the weighted average of the IR ratios for each gene
weighted_avg_ir_ratios = np.average(ir_data, axis=1, weights=weights)

# Normalize the weighted average scores to a range between 0 and 1
min_score = np.min(weighted_avg_ir_ratios)
max_score = np.max(weighted_avg_ir_ratios)
normalized_scores = (weighted_avg_ir_ratios - min_score) / (max_score - min_score)

# Save the scores as a pandas DataFrame
scores_df = pd.DataFrame(normalized_scores, columns=['score'], index=ir_data.index)

# Define the threshold for filtering
ir_threshold = 0.1

# Calculate the proportion of samples above the threshold for each gene
proportion_above_threshold = (ir_data >= ir_threshold).sum(axis=1) / len(ir_data.columns)

# Define the proportion threshold for filtering
proportion_threshold = 0.8

# Filter the genes based on the proportion_threshold
filtered_ir_data = ir_data[proportion_above_threshold >= proportion_threshold]

# Add the score column to the filtered_ir_data
filtered_ir_data_with_scores = pd.concat([filtered_ir_data, scores_df.loc[filtered_ir_data.index]], axis=1)

# Combine the metadata_columns with the filtered_ir_data_with_scores
final_filtered_data = pd.concat([metadata_columns.loc[filtered_ir_data_with_scores.index], filtered_ir_data_with_scores], axis=1)

first_four_columns = depth_df_kfo.iloc[:, :4]
remaining_columns = depth_df_kfo.iloc[:, 4:-2]

# Apply CPM normalization to each sample (column) in the dataframe
normalized_depths_kfo = remaining_columns.apply(calculate_cpm, axis=0)
medians = normalized_depths_kfo.median(axis=1)
normalized_depth_total = pd.concat([first_four_columns, normalized_depths_kfo],axis=1)

normalized_depth_total['medians'] = medians
normalized_depth_total = normalized_depth_total[normalized_depth_total['medians']!=0]


##### Scatter-Plot #####
# Assign colors and transparencies based on the gene presence in the lists and their scores
colors = []
transparencies = []
green_list = []

other_colors = []
other_transparencies = []
data_labels = []

for _, row in kfo_main_graph.iterrows():
    gene = row['Name']
    start = row['Start']
    end = row['End']
    
    if any((data_export_final_irs['Name'] == gene) & (data_export_final_irs['Start'] == start) & (data_export_final_irs['End'] == end)):
        score = data_export_final_irs.loc[(data_export_final_irs['Name'] == gene) & (data_export_final_irs['Start'] == start) & (data_export_final_irs['End'] == end), 'score'].iloc[0]
        if gene.split('/')[0] in gene_list_muts:
            green_list.append(gene.split('/')[0])

            other_colors.append('green')
            other_transparencies.append(1)
        else:
            other_colors.append('blue')
            other_transparencies.append(1)
    elif gene.split('/')[0] in gene_list_muts:
        other_colors.append('purple')
        other_transparencies.append(1)
    else:
        other_colors.append('gray')
        other_transparencies.append(1)
    
    if other_colors and other_colors[-1] == 'green':
        data_labels.append(gene.split('/')[0])
    else:
        data_labels.append(None)


# Convert the colors and transparencies to RGBA format
rgba_colors = np.zeros((len(kfo_main_graph), 4))

for i, (color, transparency) in enumerate(zip(other_colors, other_transparencies)):
    if color == 'blue':
        rgba_colors[i] = (0, 0, 1, transparency)
    elif color == 'green':
        rgba_colors[i] = (0, 1, 0, transparency)
    elif color == 'purple':
        rgba_colors[i] = (0.5, 0, 0.5, transparency)
    else:  # gray
        rgba_colors[i] = (0.5, 0.5, 0.5, transparency)

# Create a scatter plot with colored dots, increased resolution, and smaller dot sizes
plt.figure(figsize=(10, 6), dpi=250)

gray_indices = [i for i, color in enumerate(colors) if color == 'gray']
other_indices = [i for i, color in enumerate(colors) if color != 'gray']

x_values = kfo_main_graph['medians_x']
y_values = np.log2(kfo_main_graph['medians_y'])

gray_mask = (rgba_colors == [0.5, 0.5, 0.5, 1.0]).all(axis=1)
gray_points = kfo_main_graph[gray_mask]
x_gray = gray_points['medians_x']
y_gray = np.log2(gray_points['medians_y'])
plt.scatter(x_gray, y_gray, c=rgba_colors[gray_mask], s=5)

# Boolean mask for other points
other_mask = ~gray_mask
other_points = kfo_main_graph[other_mask]
x_other = other_points['medians_x']
y_other = np.log2(other_points['medians_y'])
plt.scatter(x_other, y_other, c=rgba_colors[other_mask], s=5)


rgba_other_colors = rgba_colors[other_mask]


for idx, (i, row) in enumerate(other_points.iterrows()):
    gene = row['Name'].split('/')[0]
    x = row['medians_x']
    y = np.log2(row['medians_y'])
    color = rgba_other_colors[idx]

    if color[0] == 0 and color[1] == 1 and color[2] == 0:  # Check if the color is purple
        point_size = 20  # Increase the size of purple points
    else:
        point_size = 5


# Draw purple points last to put them in front of all other points
for idx, (gene, color, x, y) in enumerate(zip(other_points['Name'].str.split('/').str[0], rgba_other_colors, other_points['medians_x'], np.log2(other_points['medians_y']))):
    if color[0] == 0 and color[1] == 1 and color[2] == 0:  # Check if the color is purple
        plt.scatter(x, y, c=[color], s=20)


# Add a legend with color labels
blue_patch = plt.Line2D([], [], color='blue', marker='o', markersize=5, label='IR Events Shared Public PDAC and CRU-5002', linestyle='')
green_patch = plt.Line2D([], [], color='purple', marker='o', markersize=5, label='High Intronic Mutation Events (Total-Prox-SS)', linestyle='')
purple_patch = plt.Line2D([], [], color='green', marker='o', markersize=5, label='Mutated and Shared IR Genes', linestyle='')
gray_patch = plt.Line2D([], [], color='gray', marker='o', markersize=5, label='IRs in CRU-5002', linestyle='')
plt.legend(handles=[blue_patch, green_patch, purple_patch, gray_patch])

plt.xlabel('Median IR Values')
plt.ylabel('Median log2(CPM-normalized Read Depths)')
plt.title('Scatter Plot of Median CPM-normalized Read Depths vs. Median IR Values (Total)')
plt.show()


