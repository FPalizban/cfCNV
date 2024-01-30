# cfCNV
After working with the data and understanding the data type, we have decided to apply a benchmarking analysis to identify the optimum pipeline for CNV calling in the area of long read cfDNA data. To fulfill this objective we need to find a proper dataset along with most related CNV callers. 
Here is the table for the datasets:



Here is the table of the selected tools:


Copy number variation (CNV) calling in the context of long-read circulating cell-free DNA (cfDNA) data involves identifying alterations in the number of copies of genomic segments in the cfDNA samples. Long-read sequencing technologies, such as those provided by platforms like Oxford Nanopore or PacBio, offer advantages in resolving complex genomic regions, including those with repetitive elements. Here's a general guide on CNV calling in long-read cfDNA data:

Preprocessing:
Read Alignment:

Align long reads to a reference genome using long-read aligners (e.g., Minimap2 for Oxford Nanopore reads or NGMLR for PacBio reads).
Quality control and filtering are crucial to ensure only high-quality alignments are considered.
Data Normalization:

Normalize read counts across samples to account for differences in sequencing depth.
CNV Calling:
Depth-of-Coverage Analysis:

Utilize tools that estimate the depth of coverage at each genomic position.
Deviations from the expected coverage can indicate CNVs.
Read Depth-Based CNV Detection:

Algorithms like CNVnator or FREEC can be adapted for long-read data.
These tools analyze changes in read depth across the genome to identify CNVs.
Split-Read Analysis:

Long reads can span the breakpoints of CNVs.
Tools like Sniffles (for structural variant detection) can be useful.
Assembly-Based Approaches:

Perform de novo assembly of long reads to identify structural variations, including CNVs.
Tools like Canu or Flye can be used for assembly.
Combine Multiple Signals:

Integrate multiple signals, including read depth, split-reads, and assembly information, to increase sensitivity and accuracy.
Post-Processing:
Filtering:

Apply filters to reduce false positives.
Consider removing CNVs in regions prone to mapping errors or artifacts.
Annotation:

Annotate detected CNVs with genomic features to prioritize biologically relevant events.
Validation:
Validation Experiments:

Validate predicted CNVs using orthogonal methods, such as quantitative PCR (qPCR) or digital droplet PCR (ddPCR).
Comparison with Short-Read Data:

Validate CNVs using data generated from short-read sequencing platforms.
Challenges and Considerations:
Noise and Errors:

Long-read data may have higher error rates. Quality filtering is crucial.
Computational Resources:

Analyzing long-read data requires significant computational resources.
Integration with Other Omics Data:

Consider integrating CNV data with other omics data for a comprehensive analysis.
Biological Interpretation:

Interpret detected CNVs in the context of the specific biological question or disease under investigation.
