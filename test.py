import pysam
from random import randint, random

def overlap(region1, region2):
    """
    Check if two regions overlap.
    region1=(CHR, start, end)
    """
    chr1 = region1[0] if region1[0].startswith("chr") else f"chr{region1[0]}"
    chr2 = region2[0] if region2[0].startswith("chr") else f"chr{region2[0]}"

    if chr1 == chr2:  # Check if chromosomes are the same
        # Check if regions overlap based on start and end positions
        return (region1[1] >= region2[1] and region1[1] <= region2[2]) or \
               (region2[1] >= region1[1] and region2[1] <= region1[2])
    else:
        return False  # Different chromosomes, no overlap

def remove_overlapping_regions(random_regions, regions_to_exclude):
    """
    Remove regions from random_regions that overlap with any region in regions_to_exclude.
    """
    non_overlapping_regions = []

    for random_region in random_regions:
        # Check if the current region overlaps with any region in regions_to_exclude
        overlapping = False
        for exclude_region in regions_to_exclude:
            if overlap(random_region, exclude_region):
                overlapping = True
                break

        # Add the region to non_overlapping_regions if it does not overlap
        if not overlapping:
            non_overlapping_regions.append(random_region)

    return non_overlapping_regions

def simulate_cnv(bam_input, bam_output, num_regions, output_regions_file, regions_to_exclude=[]):
    # Open input BAM file
    input_bam = pysam.AlignmentFile(bam_input, 'rb')

    # Open output BAM file
    output_bam = pysam.AlignmentFile(bam_output, 'wb', header=input_bam.header)

    # List to store randomly generated regions
    random_regions = []

    # Iterate through specified number of random regions
    for _ in range(num_regions):
        while True:
            # Choose a random chromosome from 1 to 22
            chrom = str(randint(1, 22))

            # Check if the randomly generated chromosome is present in the BAM file references
            if chrom not in input_bam.references:
                # Handle the case where the chromosome is not present
                print(f"Chromosome {chrom} not found in BAM file references. Trying another random chromosome.")
                continue

            # Get the index of the chromosome in the references
            chrom_index = input_bam.references.index(chrom)

            # Check if there are any reads mapped to the chosen chromosome
            if input_bam.count(chrom) == 0:
                print(f"No reads mapped to chromosome {chrom}. Trying another random chromosome.")
                continue

            start = randint(1, input_bam.lengths[chrom_index])  # Random start position
            end = min(start + randint(1000000, 5000000), input_bam.lengths[chrom_index])  # Random end position

            # Simulate copy number variation as duplication or homozygous deletion
            if random() < 0.5:
                # Simulate homozygous deletion
                copy_number = 0
            else:
                # Simulate copy number variation as duplication (values 2 to 5)
                copy_number = randint(2, 5)

            # Append the generated region to the list
            random_regions.append((chrom, start, end, copy_number))

            # Remove all reads mapped to the specified region if copy_number is 0
            if copy_number == 0:
                # Skip all reads in the specified region
                continue
            else:
                # Initialize index for appending an increasing number to the read names
                idx = 1
                # Iterate through reads and modify read positions in the specified region
                for read in input_bam.fetch(chrom, start, end):
                    new_read_name = f"{read.query_name}-{idx}"  # Append an increasing number
                    idx += 1  # Increment the index for the next read
                    read.query_name = new_read_name
                    output_bam.write(read)

            break  # Break the while loop if a valid region is generated

    # Remove overlapping regions with regions_to_exclude
    random_regions = remove_overlapping_regions(random_regions, regions_to_exclude)

    # Sort the final sample by chromosome coordinates
    random_regions.sort(key=lambda x: (x[0], x[1]))

    # Write the list of non-overlapping randomly generated regions to a file
    with open(output_regions_file, 'w') as regions_file:
        for region in random_regions:
            regions_file.write(f"{region}\n")

    # Write the remaining reads outside the simulated regions
    for chrom in input_bam.references:
        if chrom.isdigit() and 1 <= int(chrom) <= 22:
            # Check if it is a chromosome from 1 to 22
            if (chrom, 1, 1, 1) not in random_regions:  # we don't duplicate the reads from simulated regions
                for read in input_bam.fetch(chrom):
                    output_bam.write(read)

    # Close BAM files
    input_bam.close()
    output_bam.close()

    return random_regions

# Example usage with configurable parameters
input_bam_file = "/Users/apple/Documents/unil/simulation/IH01_chr22_no-black.bam"
output_bam_file = "/Users/apple/Documents/unil/simulation/IH01_chr22_no-black-simulated1.bam"
num_random_regions = 3
output_regions_file = "/Users/apple/Documents/unil/simulation/random_regions1-IH01_chr22_no-black.txt"
regions_to_exclude_file = "/Users/apple/Documents/unil/simulation/hg19-blacklist.v2.txt"

# Read regions to exclude from a file (format: chrom start end copy_number)
regions_to_exclude = []
with open(regions_to_exclude_file, 'r') as exclude_file:
    for line in exclude_file:
        parts = line.strip().split()
        if len(parts) == 4:  # Ensure the line has the expected number of elements
            regions_to_exclude.append((parts[0], int(parts[1]), int(parts[2]), int(parts[3])))
        else:
            print(f"Skipping invalid line in {regions_to_exclude_file}: {line}")

random_regions_list = simulate_cnv(input_bam_file, output_bam_file, num_random_regions, output_regions_file, regions_to_exclude)

# Print the list of randomly generated non-overlapping regions
print("Randomly generated non-overlapping regions:")
for region in random_regions_list:
    print(region)
