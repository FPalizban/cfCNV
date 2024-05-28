import csv

def read_regions_from_file1(filename):
    """
    Read region information from a text file with a specified format.
    :param filename: Name of the text file containing region information
    :return: List of tuples (chromosome, start, end, index) representing regions
    """
    regions = []
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        for row in csv_reader:
            if len(row) >= 3:
                chromosome, start, end = row[:3]
                regions.append((chromosome, int(start), int(end)))
    return regions

def read_regions_from_file2(filename):
    """
    Read region information from a text file with a specified format.
    :param filename: Name of the text file containing region information
    :return: List of tuples (chromosome, start, end) representing regions
    """
    regions = []
    with open(filename, 'r') as file:
        for line in file:
            parts = eval(line.strip())
            if len(parts) >= 3:
                chromosome, start, end = parts[:3]
                regions.append((chromosome, int(start), int(end)))
    return regions

def calculate_overlap(region1, region2):
    chromosome1, start1, end1= region1
    chromosome2, start2, end2= region2
    # Check if chromosomes match
    if chromosome1 != chromosome2:
        return 0.0

    # Calculate the intersection and union of the coordinates
    intersection = max(0, min(end1, end2) - max(start1, start2))
    union = max(end1, end2) - min(start1, start2)

    overlap_ratio = intersection / (end2 - start2) if union > 0 else 0
    return overlap_ratio

def find_overlapping_regions(region_list, region_true, threshold=0.5):
    overlapping_indices = []

    for i, region1 in enumerate(region_list):
        for j, region2 in enumerate(region_true):
            overlap_ratio = calculate_overlap(region1, region2)

            # Check if the overlap is greater than or equal to the threshold
            if overlap_ratio >= threshold:
                overlapping_indices.append((j, i))

    return overlapping_indices

# Example usage:
file1 = '/Users/apple/Documents/unil/results/simulations/IH01-sim-s1-qdna.txt'
file2 = '/Users/apple/Documents/unil/results/simulations/random_regions-test_IH01.txt'

regions1 = read_regions_from_file1(file1)
regions_true = read_regions_from_file2(file2)

overlapping_indices = find_overlapping_regions(regions1, regions_true, threshold=0.5)



unique_values = set(item[0] for item in overlapping_indices)
count_unique_values = len(unique_values)

TP= count_unique_values
FP= len(regions1) - len(set(item[1] for item in overlapping_indices))

recall = TP / len(regions_true)
precision = TP / (TP + FP)

#######################################################################

####### bar plot for each sample ########
import matplotlib.pyplot as plt

# Plotting
labels = ['True Positives', 'False Positives']
values = [TP, FP]

plt.bar(labels, values, color=['green', 'red'])
plt.title('True Positives and False Positives of spectree1000000 with 50%')
plt.xlabel('Categories')
plt.ylabel('Count')
plt.show()


###### roc cutve for each bin size in each tool ########
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve

# Assuming you have your TP and FP data in two lists
TP = [9, 4, 7, 5]
FP = [0, 0, 0, 0]

# Calculate the true positive rate (TPR) and false positive rate (FPR)
TPR = [t / max(TP) for t in TP]
FPR = [f / max(FP) for f in FP]

# Plot the ROC curve
plt.plot(TP, FP, marker='o')
plt.title('ROC Curve of Spectre-1000000 with 50%')
plt.xlabel('False Positive')
plt.ylabel('True Positive ')
plt.grid(True)
plt.show()



####### multiple roc curve for each bin size with multiple tools ########

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve

# Assuming you have TP and FP data for multiple tools
TP_tool1 = [0.2, 0.4, 0.6, 0.7, 0.8]
FP_tool1 = [0.1, 0.2, 0.3, 0.5, 0.6]

TP_tool2 = [0.1, 0.3, 0.5, 0.6, 0.7]
FP_tool2 = [0.2, 0.3, 0.4, 0.5, 0.7]

# Calculate TPR and FPR for each tool
TPR_tool1 = [t / max(TP_tool1) for t in TP_tool1]
FPR_tool1 = [f / max(FP_tool1) for f in FP_tool1]

TPR_tool2 = [t / max(TP_tool2) for t in TP_tool2]
FPR_tool2 = [f / max(FP_tool2) for f in FP_tool2]

# Plot the ROC curves for each tool
plt.plot(FPR_tool1, TPR_tool1, marker='o', label='Tool 1')
plt.plot(FPR_tool2, TPR_tool2, marker='o', label='Tool 2')

plt.title('ROC Curve')
plt.xlabel('False Positive Rate (FPR)')
plt.ylabel('True Positive Rate (TPR)')
plt.grid(True)
plt.legend()
plt.show()


############# precision and recall for each tool with different bin sizes #########

import matplotlib.pyplot as plt

# Example precision and recall values
experiments = ['Spectre-1000', 'Spectre-100000', 'Spectre-1000000']
precision_values = [0.8701298701298701, 0.8604651162790697, 1.0]
recall_values = [0.34986945169712796, 0.38642297650130547, 0.02349869451697128]

plt.plot(experiments, precision_values, marker='o', label='Precision')
plt.plot(experiments, recall_values, marker='o', label='Recall')

plt.title('Precision and Recall for Spectre with different bin sizes with 50% overlap')
plt.xlabel('Bin sizes')
plt.ylabel('Score')
plt.legend()
plt.grid(True)
plt.show()


import matplotlib.pyplot as plt

# Example precision and recall values
experiments = ['Spectre-1000', 'Spectre-100000', 'Spectre-1000000']
precision_values = [0.8701298701298701, 0.8604651162790697, 1.0]
recall_values = [0.34986945169712796, 0.38642297650130547, 0.02349869451697128]

plt.figure(figsize=(8, 6))
plt.plot(recall_values, precision_values, marker='o', linestyle='-', color='b', label='Precision-Recall Curve')

# Annotate points with experiment names
for i, txt in enumerate(experiments):
    plt.annotate(txt, (recall_values[i], precision_values[i]), textcoords="offset points", xytext=(10,-5), ha='center', fontsize=10)

# Highlight points with larger markers
plt.scatter(recall_values, precision_values, color='r', s=100, label='Bin size')

plt.title('Precision-Recall Curve for Spectre with Different Bin Sizes')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


#######################################################################
import matplotlib.pyplot as plt

# Example precision and recall values with replicates
experiments = ['Spectre-1000', 'Spectre-100000', 'Spectre-1000000']
precision_values = [[0.8701298701298701,0.8404907975460123, 0.88268156424581, 0.8819875776397516], [0.8604651162790697,0.90625,0.8895027624309392,0.8947368421052632], [1, 1,1 ,1.0]]
recall_values = [[0.34986945169712796,0.3366093366093366, 0.3691588785046729, 0.36883116883116884], [0.38642297650130547,0.4275184275184275,0.37616822429906543,0.3974025974025974], [0.02349869451697128, 0.009828009828009828, 0.016355140186915886,0.012987012987012988]]

plt.figure(figsize=(8, 6))

# Plot precision-recall curves for each experiment
for i, experiment in enumerate(experiments):
    plt.plot(recall_values[i], precision_values[i], marker='o', linestyle='-', label=experiment)

plt.title('Precision-Recall Curve for Spectre with Different Bin Sizes and Replicates')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# bar plot
import numpy as np
import matplotlib.pyplot as plt

# Example precision and recall values with replicates
experiments = ['Spectre-1000', 'Spectre-100000', 'Spectre-1000000']
precision_values = [
    [0.8701298701298701, 0.8404907975460123, 0.88268156424581, 0.8819875776397516],
    [0.8604651162790697, 0.90625, 0.8895027624309392, 0.8947368421052632],
    [1, 1, 1, 1.0]
]
recall_values = [
    [0.34986945169712796, 0.3366093366093366, 0.3691588785046729, 0.36883116883116884],
    [0.38642297650130547, 0.4275184275184275, 0.37616822429906543, 0.3974025974025974],
    [0.02349869451697128, 0.009828009828009828, 0.016355140186915886, 0.012987012987012988]
]

# Calculate mean precision and recall values along with their standard deviations
mean_precision_values = [np.mean(vals) for vals in precision_values]
std_precision_values = [np.std(vals) for vals in precision_values]

mean_recall_values = [np.mean(vals) for vals in recall_values]
std_recall_values = [np.std(vals) for vals in recall_values]

# Plot bar plots with error bars
plt.figure(figsize=(10, 6))

bar_width = 0.35
index = np.arange(len(experiments))

plt.bar(index, mean_precision_values, bar_width, label='Precision', yerr=std_precision_values)
plt.bar(index + bar_width, mean_recall_values, bar_width, label='Recall', yerr=std_recall_values)

plt.xlabel('Experiments')
plt.ylabel('Scores')
plt.title('Mean Precision and Recall with Standard Deviation')
plt.xticks(index + bar_width / 2, experiments)
plt.legend()

plt.tight_layout()
plt.show()


#########################################
experiments = ['ichor-40-1000', 'ichor-40-100000', 'ichor-40-1000000']

precision_values = [
    [0.7701298701298701, 0.6404907975460123, 0.78268156424581, 0.7819875776397516],
    [0.6604651162790697,0.8599221789883269,0.6885245901639344, 0.737012987012987],
    [0.67, 0.80, 0.70, 0.80]
]
recall_values = [
    [0.34986945169712796, 0.2366093366093366, 0.4691588785046729, 0.36883116883116884],
    [0.38642297650130547, 0.5163551401869159, 0.5454545454545454, 0.535377358490566],
    [0.03349869451697128, 0.07828009828009828, 0.036355140186915886, 0.0312987012987012988]
]

mean_precision = []
mean_recall = []

for i in range(len(experiments)):
    mean_precision.append(sum(precision_values[i]) / len(precision_values[i]))
    mean_recall.append(sum(recall_values[i]) / len(recall_values[i]))

print("Experiment\tMean Precision\tMean Recall")
for i in range(len(experiments)):
    print(f"{experiments[i]}\t{mean_precision[i]:.6f}\t\t{mean_recall[i]:.6f}")



# Calculate mean precision and recall values along with their standard deviations
mean_precision_values = [np.mean(vals) for vals in precision_values]
std_precision_values = [np.std(vals) for vals in precision_values]

mean_recall_values = [np.mean(vals) for vals in recall_values]
std_recall_values = [np.std(vals) for vals in recall_values]

# Plot bar plots with error bars
plt.figure(figsize=(10, 6))

bar_width = 0.35
index = np.arange(len(experiments))

plt.bar(index, mean_precision_values, bar_width, label='Precision', yerr=std_precision_values)
plt.bar(index + bar_width, mean_recall_values, bar_width, label='Recall', yerr=std_recall_values)

plt.xlabel('Experiments')
plt.ylabel('Scores')
plt.title('Mean Precision and Recall with Standard Deviation')
plt.xticks(index + bar_width / 2, experiments)
plt.legend()

plt.tight_layout()
plt.show()



#######################################
#########################################
experiments = ['QDNAseq100000', 'QDNAseq-1000000']

precision_values = [
    [0.6604651162790697,0.8599221789883269,0.6885245901639344, 0.737012987012987],
    [1.0, 1.0,1.0, 1.0]
]
recall_values = [
    [0.38642297650130547, 0.5163551401869159, 0.5454545454545454, 0.535377358490566],
    [0.1957547169811321, 0.2805194805194805, 0.3037383177570093, 0.25307125307125306]
]

mean_precision = []
mean_recall = []

for i in range(len(experiments)):
    mean_precision.append(sum(precision_values[i]) / len(precision_values[i]))
    mean_recall.append(sum(recall_values[i]) / len(recall_values[i]))

print("Experiment\tMean Precision\tMean Recall")
for i in range(len(experiments)):
    print(f"{experiments[i]}\t{mean_precision[i]:.6f}\t\t{mean_recall[i]:.6f}")



# Calculate mean precision and recall values along with their standard deviations
mean_precision_values = [np.mean(vals) for vals in precision_values]
std_precision_values = [np.std(vals) for vals in precision_values]

mean_recall_values = [np.mean(vals) for vals in recall_values]
std_recall_values = [np.std(vals) for vals in recall_values]

# Plot bar plots with error bars
plt.figure(figsize=(10, 6))

bar_width = 0.35
index = np.arange(len(experiments))

plt.bar(index, mean_precision_values, bar_width, label='Precision', yerr=std_precision_values)
plt.bar(index + bar_width, mean_recall_values, bar_width, label='Recall', yerr=std_recall_values)

plt.xlabel('Experiments')
plt.ylabel('Scores')
plt.title('Mean Precision and Recall with Standard Deviation')
plt.xticks(index + bar_width / 2, experiments)
plt.legend()

plt.tight_layout()
plt.show()



###############################################
#########################################
experiments = ['cnvkit']

precision_values = [
    [0.018902162407379405, 0.02077447232840286,0.017636959575464337, 0.02105101413644745]
]
recall_values = [
    [0.3071253071253071, 0.29205607476635514, 0.2935064935064935, 0.3231132075471698]
]

mean_precision = []
mean_recall = []

for i in range(len(experiments)):
    mean_precision.append(sum(precision_values[i]) / len(precision_values[i]))
    mean_recall.append(sum(recall_values[i]) / len(recall_values[i]))

print("Experiment\tMean Precision\tMean Recall")
for i in range(len(experiments)):
    print(f"{experiments[i]}\t{mean_precision[i]:.6f}\t\t{mean_recall[i]:.6f}")



# Calculate mean precision and recall values along with their standard deviations
mean_precision_values = [np.mean(vals) for vals in precision_values]
std_precision_values = [np.std(vals) for vals in precision_values]

mean_recall_values = [np.mean(vals) for vals in recall_values]
std_recall_values = [np.std(vals) for vals in recall_values]

# Plot bar plots with error bars
plt.figure(figsize=(10, 6))

bar_width = 0.35
index = np.arange(len(experiments))

plt.bar(index, mean_precision_values, bar_width, label='Precision', yerr=std_precision_values)
plt.bar(index + bar_width, mean_recall_values, bar_width, label='Recall', yerr=std_recall_values)

plt.xlabel('Experiments')
plt.ylabel('Scores')
plt.title('Mean Precision and Recall with Standard Deviation')
plt.xticks(index + bar_width / 2, experiments)
plt.legend()

plt.tight_layout()
plt.show()



