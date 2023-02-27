import argparse
from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from ete3 import Tree
import math
import os
import random
import sys
from tqdm import tqdm

def get_options():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description='Simulate bacterial genomes from a reference and Newick tree.')
    parser.add_argument('--tree', dest='treeFile',
                        help='Path to a tree in Newick format.', required=True)
    parser.add_argument('--reference', dest='referenceFile',
                        help='Path to a bacterial reference sequence (only a single contig pls).', required=True)
    parser.add_argument('--output', dest='outputDir',
                        help='Output directory for the simulated sequences and list of SNPs.', required=True)
    parser.add_argument('--mutation-rate', dest='alpha',
                        help='The mutation rate in SNPs/site/year. Default is 1.57×10^−6 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3648787/#:~:text=The%20rate%20at%20which%20base,and%20much%20higher%20than%20that)',
                        type=float, default=0.00000157, required=False)
    parser.add_argument('--nj-tree', dest='nj_tree', action='store_true', default=False,
                        help='Write an NJ tree (Warning: it takes ages)')
    args = parser.parse_args()
    return args

def load_newick_tree(treeFile):
    """ load and return a newick tree using the ete3 module """
    # load the tree
    tree = Tree(treeFile, format=1)
    return tree

def load_reference_genome(referenceFile):
    """ load and return a reference genome in FASTA format using Bio.SeqIO and convert to a dictionary """
    # load the reference genome
    referenceGenome = SeqIO.to_dict(SeqIO.parse(open(referenceFile),'fasta'))
    # ensure the reference genome only contains 1 contig
    assert len(referenceGenome) == 1, "The reference cannot contain multiple contigs!"
    return referenceGenome

def get_snp_counts(tree,
                referenceGenome,
                alpha):
    """ returns a dictionary of the number of SNPs between each parent and offspring in the tree based on branch lengths """
    # store the number of SNPs in a dict
    snpCounts = {}
    # keep track of the mapping of topology id to node name
    nodeMapping = {}
    # get the reference genome length
    for contig in referenceGenome:
        referenceLength = len(str(referenceGenome[contig].seq))
        break
    assert referenceLength, "Could not get the length of this reference genome"
    # traverse the nodes in the graph
    for node in tree.traverse():
        # if this node is not the root node
        if not node.is_root():
            # get the time between parent and offspring in years
            evolutionaryTime = node.dist
            # calculate the number of SNPs based on the mutation rate
            numberOfSNPs =  alpha * referenceLength * evolutionaryTime
            # round the number of SNPs up to the nearest integer
            numberOfSNPs = int(math.ceil(numberOfSNPs))
            # check that this parent offspring pair has not been seen before
            assert not (node.up.get_topology_id(), node.get_topology_id()) in snpCounts, "This parent offspring pair has been seen before"
            # add the node and the number of SNPs to it's parent to the dictionary as a tuple in the format (parent, offspring)
            snpCounts[(node.up.get_topology_id(), node.get_topology_id())] = numberOfSNPs
            # if the node has a name (is a tip), then keep track of it's topology ID
            if not node.name == "":
                nodeMapping[node.get_topology_id()] = node.name
        else:
            # get the topology ID of the root node
            referenceTopologyID = node.get_topology_id()
            # store the topology ID mapping
            nodeMapping[node.get_topology_id()] = "0"
    return snpCounts, nodeMapping, referenceLength, referenceTopologyID

def get_SNP_positions(snpCounts,
                    referenceLength):
    """ return a dictionary listing the 0-based base position where we will be adding a SNP """
    # see the random number generator for consistency across runs
    random.seed(2023)
    # store the SNP positions
    SNPLocations = {}
    # iterate through each parent offspring combination
    for parentOffspringPair in snpCounts:
        # get the number of SNPs
        numberOfSnps = snpCounts[parentOffspringPair]
        # use the random module to decide the location of the SNPs
        SNPsites = [random.randint(0, referenceLength - 1) for i in range(numberOfSnps)]
        # make sure every SNP is occurring in a different position for this sample
        while not len(SNPsites) == len(set(SNPsites)):
            # use the random module to decide the location of the SNPs
            SNPsites = [random.randint(0, referenceLength - 1) for i in range(numberOfSnps)]
        # check that the number of SNPs equals the number of SNP sites
        assert numberOfSnps == len(SNPsites), "Incorrect number of SNP sites"
        # add the SNP sites to the SNPLocations dict
        SNPLocations[parentOffspringPair] = SNPsites
    return SNPLocations

def initialise_simulation_dict(referenceGenome,
                            referenceTopologyID):
    """ returns a dictionary containing the topology ID and sequence of the reference """
    # get the reference sequence without a header
    for contig in referenceGenome:
        referenceSequence = str(referenceGenome[contig].seq)
        break
    # initialise the dictionary with the reference genome information
    initialisedSimulation = {referenceTopologyID: referenceSequence}
    return initialisedSimulation

def simulate_from_reference(simulatedGenomes,
                        SNPLocations,
                        referenceTopologyID):
    """ returns a dictionary of simulated sequences across a tree """
    # keep track of the mutations that have occurred between each parent and offspring, initialising with the reference
    mutationTracking = {referenceTopologyID: []}
    # while we have not simulated all sequences
    while not len(simulatedGenomes) == len(SNPLocations) + 1:
        # iterate through the parent offspring pairs
        for pair in SNPLocations:
            # get the parent and child topology ID
            parentID = pair[0]
            childID = pair[1]
            # if we have not already simulated the child sequence and we have the parent sequence
            if parentID in simulatedGenomes:
                if not childID in simulatedGenomes:
                    # add the child offspring to the mutation tracking dict and the mutations that have occurred so far
                    mutationTracking[childID] = mutationTracking[parentID].copy()
                    # get the parent sequence
                    parentSequence = simulatedGenomes[parentID]
                    # initialise the child sequence list with the parent sequence
                    childSequence = list(parentSequence)
                    # iterate through the SNP positions
                    for pos in SNPLocations[pair]:
                        # see what the current base is in this position
                        currentBase = childSequence[pos]
                        # randomly decide which base we will mutate to
                        newBase = random.choice(["A", "C", "T", "G"])
                        # resample if the new base is the same as the previous base
                        while newBase == currentBase:
                            newBase = random.choice(["A", "C", "T", "G"])
                        # replace the base with the SNP
                        childSequence[pos] = newBase
                        # add the mutation to the mutation tracking dictionary
                        mutationTracking[childID].append({pos: (currentBase, newBase)})
                    # sort the SNP list by position
                    mutationTracking[childID] = list(sorted(mutationTracking[childID], key=lambda x: list(x.keys())[0]))
                    # convert the child sequence back to a string
                    childSequence =  "".join(childSequence)
                    # ensure the parent and child sequences are not the same
                    assert not parentSequence == childSequence, "Parent and child simulated sequences are the same"
                    assert len(parentSequence) == len(childSequence), "Off by one error"
                    # add the newly mutated child sequence to the dictionary of simulated sequences
                    simulatedGenomes[childID] = childSequence
    return simulatedGenomes, mutationTracking

def write_simulated_genomes(simulatedGenomes,
                            nodeMapping,
                            outputDir):
    """ writes a single FASTA file of the simulated sequences """
    # write out the simulated genomes
    tipGenomes = []
    # iterate through all simulated genomes
    for node in simulatedGenomes:
        # only get the sequence if the sample is a tip
        if node in nodeMapping:
            # get the sample name
            sampleName = nodeMapping[node]
            # get the sample sequence
            sampleSequence = simulatedGenomes[node]
            # add the sample name and sequence to a list
            tipGenomes.append(">" + sampleName + "\n" + sampleSequence)
    # sort the genomes by sample name
    tipGenomes = list(sorted(tipGenomes, key=lambda x: int(x.split("\n")[0].split(">")[1])))
    # write the sequences to a file
    with open(os.path.join(outputDir, "simulated_genomes.fasta"), "w") as outFasta:
        outFasta.write("\n".join(tipGenomes))
    return os.path.join(outputDir, "simulated_genomes.fasta")

def write_simulated_mutations(mutationTracking,
                            nodeMapping,
                            outputDir):
    """ writes a single txt file of the simulated mutations in each sample relative to the reference """
    tipMutations = []
    # iterate through all simulated genomes
    for node in mutationTracking:
        # only get the sequence if the sample is a tip
        if node in nodeMapping:
            # get the sample name
            sampleName = nodeMapping[node]
            # get the sample mutations
            sampleMutations = mutationTracking[node]
            # make the mutation positions more human readable
            friendlySampleMutations = [str(list(p.keys())[0]) + "\t" + "->".join(list(list(p.values())[0])) for p in sampleMutations]
            # add the sample name and mutations to a list
            tipMutations.append(">" + sampleName + "\n" + "\n".join(friendlySampleMutations))
    # sort the genomes by sample name
    tipMutations = list(sorted(tipMutations, key=lambda x: int(x.split("\n")[0].split(">")[1])))
    # write the mutations to a file
    with open(os.path.join(outputDir, "simulated_mutations.txt"), "w") as outMut:
        outMut.write("\n".join(tipMutations))

def write_nj_tree(simulatedGenomes,
                outputDir):
    """ write an NJ tree to visualise how well the initial tree structure was simulated using Bio.phylo """
    # load the genomes as an alignment as they are all the same length
    alignment = AlignIO.read(simulatedGenomes, 'fasta')
    # calculate sample distances
    sys.stderr.write("\n\tCalculating sample distances")
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    # build the NJ tree
    sys.stderr.write("\n\tConstructing NJ tree")
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.nj(distance_matrix)
    # write out the tree
    sys.stderr.write("\n\tWriting NJ tree")
    Phylo.write(tree, os.path.join(outputDir, "simulated_sample_tree.nwk"), "newick")

def main():
    # parse command line arguments
    args = get_options()
    # import the tree file
    sys.stderr.write("\nLoading tree file\n")
    tree = load_newick_tree(args.treeFile)
    # import the reference sequence
    sys.stderr.write("\nLoading reference genome\n")
    referenceGenome = load_reference_genome(args.referenceFile)
    # calculate the number of SNPs for every branch from parent to offspring
    sys.stderr.write("\nDetermining number of SNPs per node\n")
    snpCounts, nodeMapping, referenceLength, referenceTopologyID = get_snp_counts(tree,
                                                                            referenceGenome,
                                                                            args.alpha)
    # decide where we are independently adding the SNPs in each parent offspring pair
    sys.stderr.write("\nDetermining SNP positions per sample\n")
    SNPLocations = get_SNP_positions(snpCounts,
                                    referenceLength)
    # initialise a dictionary of the simulated genome with the reference
    initialisedSimulation = initialise_simulation_dict(referenceGenome,
                                                    referenceTopologyID)
    # simulate the sequence using the SNP position information
    sys.stderr.write("\nSimulating SNPs\n")
    simulatedGenomes, mutationTracking = simulate_from_reference(initialisedSimulation,
                                                                SNPLocations,
                                                                referenceTopologyID)
    # make the output directories if they doesn't exist
    if not os.path.exists(args.outputDir):
        os.mkdir(args.outputDir)
    # write out the simulated genomes to a file
    sys.stderr.write("\nWriting simulated genomes\n")
    simulatedSamplePath = write_simulated_genomes(simulatedGenomes,
                                            nodeMapping,
                                            args.outputDir)
    # write out a file listing the mutations
    sys.stderr.write("\nWriting list of mutations\n")
    write_simulated_mutations(mutationTracking,
                            nodeMapping,
                            args.outputDir)
    # write a neighborhood joining tree
    if args.nj_tree:
        sys.stderr.write("\nWriting NJ tree of simulated samples:")
        write_nj_tree(simulatedSamplePath,
                    args.outputDir)
    sys.exit(0)

if __name__ == "__main__":
    main()