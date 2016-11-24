#!/usr/bin/env python


"""
Developer: Emzo de los Santos
compare a set of clusters to MiBIG-Database using cluster distance method developed in BiG-SCAPE

Make New Data Structure for MiBIG Clusters that algorithm can use
"""
import os
from glob import glob
from functions import parsePFD,fasta_parser
from scipy.optimize import linear_sum_assignment
from Bio.SubsMat.MatrixInfo import pam250 as scoring_matrix
from Bio import pairwise2
import numpy as np
from bigscape import calculate_GK

def buildMiBIGDct(path):
    """
    Given a path to MiBIG pfs and pfd files will build a Data Structure that can be used for MiBIG comparisons, this
    includes the different domains in the MiBIG cluster, the product and the domain sequences
    :param path:
    :return: {MiBIG Cluster ID: (product,domList,domDict)}
    domDict - {pfamID:{domIdx:sequence}}
    domDict is built from getting the fasta sequences after parsePFD

    parsePFD returns:
    (clusterPfamDict,domList)

    clusterPfamDict is {Pfam ID: [(Domain ID (geneID,(start,end)) ,location)]}
    """
    # Only Build a Dictionary for Entries that have fasta files and pfdFiles

    mibigDict = {}
    pfdFiles = glob(path + '/*.pfd')
    fastaFiles = glob(path + '/*.fasta')

    pfdSet = set(x.split('.pfd')[0] for x in pfdFiles)
    fastaSet = set(x.split('.fasta')[0] for x in fastaFiles)

    clustersToProcess = pfdSet & fastaSet
    for mibigCluster in clustersToProcess:
        clusterProduct = mibigCluster.split('.')[-1]
        clusterID = mibigCluster.split(os.sep)[-1]
        mibig_pfd_handle = os.path.join(path, mibigCluster + '.pfd')
        mibig_fasta_handle = os.path.join(path,mibigCluster + '.fasta')
        clusterPfamDict,domList = parsePFD(mibig_pfd_handle)
        clusterFastaDict = fasta_parser(open(mibig_fasta_handle))
        domDict = {}
        for pfamDomain in clusterPfamDict.keys():
            domainPfamDict = {}
            for domainInstance in clusterPfamDict[pfamDomain]:
                geneID, (domStart, domEnd) = domainInstance
                domainPfamDict[domainInstance] = clusterFastaDict['>'+geneID][domStart:domEnd]
            domDict[pfamDomain] = domainPfamDict
        mibigDict[clusterID] = (clusterProduct,domList,domDict)

    return mibigDict

def cluster_mibig_distance(clusterHandle,mibigDct,anchor_domains):
    """
    calculate the distance between a specified cluster (with pfd and fasta files) to the mibig database
    :param cluster: path to cluster
    :param mibigDct: dictionary for miBIG with the correct structure (can use buildMiBIGDct to generate)
    :param anchor_domains: text file specifying anchor domains
    :return: distance_dictionary for cluster
     {mibigID:(Distance, Jaccard, DDS, GK, DDS_non_anchor, DDS_anchor, S, S_anchor)}
    """

    ## Parameters for Distance Calculation #######

    Jaccardw = 0.2
    DDSw = 0.75
    GKw = 0.05

    domain_difference_anchor, S_anchor = 0, 0
    domain_difference, S = 0, 0
    gap_open = -15
    gap_extend = -6.67
    anchorweight = 2
    nbhood = 4
    ###########################

    if os.path.isfile(os.path.join(clusterHandle + '.pfd')):
        cluster_pfd_handle = os.path.join(clusterHandle + '.pfd')
    else:
        print "No PFD file for %s" % clusterHandle
        raise Exception
    if os.path.isfile(os.path.join(clusterHandle + '.fasta')):
        cluster_fasta_handle = os.path.join(clusterHandle + '.fasta')
    else:
        print "No Fasta file for %s" % clusterHandle
        raise Exception

    cluster_fasta_dict = fasta_parser(open(cluster_fasta_handle))
    cluster, clusterDomList = parsePFD(cluster_pfd_handle)

    clusterDomSet = set(clusterDomList)

    cluster_dist_dict = {}

    for (mibigCluster,(product,mibigDomList,mibigDomDict)) in mibigDct.iteritems():
        print mibigCluster
        mibigDomSet = set(mibigDomDict.keys())

        intersect = clusterDomSet.intersection(mibigDomSet)
        not_intersect = clusterDomSet.symmetric_difference(mibigDomSet)

        for unshared_domain in not_intersect:  # no need to look at seq identity, since these domains are unshared
            # for each occurence of an unshared domain do domain_difference += count of domain and S += count of domain
            unshared_occurrences = []
            if unshared_domain in clusterDomSet:
                unshared_occurrences = cluster[unshared_domain]
            else:
                unshared_occurrences = mibigDomDict[unshared_domain]

            # don't look at domain version, hence the split
            if unshared_domain.split(".")[0] in anchor_domains:
                domain_difference_anchor += len(unshared_occurrences)
            else:
                domain_difference += len(unshared_occurrences)
        S = domain_difference  # can be done because it's the first use of these
        S_anchor = domain_difference_anchor

        #These are the cases for all the shared domains
        for shared_domain in intersect:
            # First we need to get the sequences and generate the domain dictionary for the cluster
            clusterDomIdx = cluster[shared_domain]
            clusterPfamDict = dict()
            for domain in clusterDomIdx:
                geneID,(domStart,domEnd) = domain
                clusterPfamDict[domain] = cluster_fasta_dict['>'+geneID][domStart:domEnd]

            mibigPfamDict = mibigDomDict[shared_domain]

            clusSize = len(clusterPfamDict)
            mibigSize = len(mibigPfamDict)

            scoreMatrix = np.ndarray((clusSize,mibigSize))
            # Get all of the instances of the domain from cluster A and cluster B and populate a similarity matrix
            for i,domA in enumerate(clusterPfamDict.keys()):
                for j,domB in enumerate(mibigPfamDict.keys()):

                    seqA = clusterPfamDict[domA]
                    seqB = mibigPfamDict[domB]

                    # this will guarantee that youll always get symmetry with comparisons by ordering it alphabetically
                    # but we'll lose accuracy and won't be able to recover the alignments

                    if seqA > seqB:
                        seqB,seqA = seqA,seqB
                    # Using BioPython
                    alignScore = pairwise2.align.globalds(seqA, seqB, scoring_matrix, gap_open, gap_extend)
                    # calculate percent similarity from best scoring alignment
                    bestAlignment = alignScore[0]
                    alignA = bestAlignment[0]
                    alignB = bestAlignment[1]
                    posCtr = 0.
                    for (a,b) in zip(alignA,alignB):
                        if a == '-' or b == '-':
                            pass
                        else:
                            if a == b:
                                posCtr += 1
                    # score is 1 - % identity since we want to define a distance to minimize
                    scoreMatrix[i,j] = 1 - posCtr/len(alignA)
            # used scipy's linear_sum_assigment to do the hungarian algorithm, haven't tested Munkres for behaviour

            pairings = [(x,y) for x,y in zip(*linear_sum_assignment(scoreMatrix)) if (x<clusSize) and (y<mibigSize)]
            # total distance from all of the paired domains
            accumulated_distance = sum(scoreMatrix[a] for a in pairings)
            # to get the total sequence distance you need to add the unpaired domains
            sum_seq_dist = accumulated_distance + abs(clusSize - mibigSize)
            # update the DDS or DDS_anchor and total domain counts depending on whether or not the domain is an anchor domain
            if shared_domain.split(".")[0] in anchor_domains:
                S_anchor += max(clusSize,mibigSize)
                domain_difference_anchor += sum_seq_dist
            else:
                S += max(clusSize,mibigSize)
                domain_difference += sum_seq_dist

        if S_anchor != 0 and S != 0:
            DDS_non_anchor = domain_difference / float(S)
            DDS_anchor = domain_difference_anchor / float(S_anchor)

            # Calculate proper, proportional weight to each kind of domain
            non_anchor_prct = S / float(S + S_anchor)
            anchor_prct = S_anchor / float(S + S_anchor)

            # boost anchor subcomponent and re-normalize
            non_anchor_weight = non_anchor_prct / (anchor_prct * anchorweight + non_anchor_prct)
            anchor_weight = anchor_prct * anchorweight / (anchor_prct * anchorweight + non_anchor_prct)

            # Use anchorweight parameter to boost percieved rDDS_anchor
            DDS = (non_anchor_weight * DDS_non_anchor) + (anchor_weight * DDS_anchor)

        elif S_anchor == 0:
            DDS_non_anchor = domain_difference / float(S)
            DDS_anchor = 0.0
            DDS = DDS_non_anchor

        else:  # only anchor domains were found
            DDS_non_anchor = 0.0
            DDS_anchor = domain_difference_anchor / float(S_anchor)
            DDS = DDS_anchor

        DDS = 1 - DDS  # transform into similarity
        # GK INDEX
        #  calculate the Goodman-Kruskal gamma index

        clusterDomListr = [item for item in clusterDomList]
        clusterDomListr.reverse()

        GK = max([calculate_GK(clusterDomList, mibigDomList, nbhood), calculate_GK(clusterDomListr, mibigDomList, nbhood)])
        Jaccard = len(intersect) / float(len(set(clusterDomList)) + len(set(mibigDomList)) - len(intersect))

        Distance = 1 - (Jaccardw * Jaccard) - (DDSw * DDS) - (GKw * GK)

        cluster_dist_dict[mibigCluster] = (Distance, Jaccard, DDS, GK, DDS_non_anchor, DDS_anchor, S, S_anchor)
    return cluster_dist_dict

def cluster_distance_pairwise_align(A, B, outputdir,  anchor_domains):
    try:
        """Compare two clusters using information on their domains, and the sequences of the domains this
        version of the script does not need the BGCs or DMS dictionary but requires the pfs and pfd files
        """
        ## Check and load required files for the cluster
        if os.path.isfile(os.path.join(outputdir, A + '.pfd')):
            clusterA_pfd_handle = os.path.join(outputdir, A + '.pfd')
        else:
            print "No PFD file for %s" % A
            raise Exception

        if os.path.isfile(os.path.join(outputdir, B + '.pfd')):
            clusterB_pfd_handle = os.path.join(outputdir, B + '.pfd')
        else:
            print "No PFD file for %s" % B
            raise Exception

        if os.path.isfile(os.path.join(outputdir, A + '.fasta')):
            clusterA_fasta_handle = os.path.join(outputdir, A + '.fasta')
            clusterA_fasta_dict = fasta_parser(open(clusterA_fasta_handle))
        else:
            print "No Fasta file for %s" % A
            raise Exception
        if os.path.isfile(os.path.join(outputdir, B + '.fasta')):
            clusterB_fasta_handle = os.path.join(outputdir, B + '.fasta')
            clusterB_fasta_dict = fasta_parser(open(clusterB_fasta_handle))
        else:
            print "No Fasta file for %s" % B
            raise Exception

        clusterA,A_list = parsePFD(clusterA_pfd_handle)
        clusterB,B_list = parsePFD(clusterB_pfd_handle)

        A_domlist = set(A_list)
        B_domlist = set(B_list)

        intersect = set(A_domlist).intersection(B_domlist)
        not_intersect = set(A_domlist).symmetric_difference(set(B_domlist))

        # JACCARD INDEX
        Jaccard = len(intersect) / float(len(set(A_domlist)) + len(set(B_domlist)) - len(intersect))

        # DDS INDEX
        # domain_difference: Difference in sequence per domain. If one cluster doesn't have a domain at all, but the other does,
        # this is a sequence difference of 1. If both clusters contain the domain once, and the sequence is the same, there is a seq diff of 0.
        # S: Max occurence of each domain
        domain_difference_anchor, S_anchor = 0, 0
        domain_difference, S = 0, 0
        gap_open = -15
        gap_extend = -6.67

        pair = ""  # pair of clusters to access their sequence identity

        # Case 1
        for unshared_domain in not_intersect:  # no need to look at seq identity, since these domains are unshared
            # for each occurence of an unshared domain do domain_difference += count of domain and S += count of domain
            unshared_occurrences = []
            if unshared_domain in A_domlist:
                unshared_occurrences = clusterA[unshared_domain]
            else:
                unshared_occurrences = clusterB[unshared_domain]

            # don't look at domain version, hence the split
            if unshared_domain.split(".")[0] in anchor_domains:
                domain_difference_anchor += len(unshared_occurrences)
            else:
                domain_difference += len(unshared_occurrences)
        S = domain_difference  # can be done because it's the first use of these
        S_anchor = domain_difference_anchor

        #These are the cases for all the shared domains
        for shared_domain in intersect:
            # First we want to index all the sequences of the domains in a dictionary
            domIdxA = clusterA[shared_domain]
            domIdxB = clusterB[shared_domain]
            domDictA = dict()
            for domain in domIdxA:
                geneID,(domStart,domEnd) = domain
                domDictA[domain] = clusterA_fasta_dict['>'+geneID][domStart:domEnd]

            domDictB = dict()
            for domain in domIdxB:
                geneID,(domStart,domEnd) = domain
                domDictB[domain] = clusterB_fasta_dict['>'+geneID][domStart:domEnd]

            clusAsize = len(domDictA)
            clusBsize = len(domDictB)

            scoreMatrix = np.ndarray((clusAsize,clusBsize))
            # Get all of the instances of the domain from cluster A and cluster B and populate a similarity matrix
            for i,domA in enumerate(domIdxA):
                for j,domB in enumerate(domIdxB):

                    seqA = domDictA[domA]
                    seqB = domDictB[domB]

                    # this will guarantee that youll always get symmetry with comparisons by ordering it alphabetically
                    # but we'll lose accuracy and won't be able to recover the alignments

                    if seqA > seqB:
                        seqB,seqA = seqA,seqB
                    # Using BioPython
                    alignScore = pairwise2.align.globalds(seqA, seqB, scoring_matrix, gap_open, gap_extend)
                    # calculate percent similarity from best scoring alignment
                    bestAlignment = alignScore[0]
                    alignA = bestAlignment[0]
                    alignB = bestAlignment[1]
                    posCtr = 0.
                    for (a,b) in zip(alignA,alignB):
                        if a == '-' or b == '-':
                            pass
                        else:
                            if a == b:
                                posCtr += 1
                    # score is 1 - % identity since we want to define a distance to minimize
                    scoreMatrix[i,j] = 1 - posCtr/len(alignA)
            # used scipy's linear_sum_assigment to do the hungarian algorithm, haven't tested Munkres for behaviour

            pairings = [(x,y) for x,y in zip(*linear_sum_assignment(scoreMatrix)) if (x<clusAsize) and (y<clusBsize)]
            # total distance from all of the paired domains
            accumulated_distance = sum(scoreMatrix[a] for a in pairings)
            # to get the total sequence distance you need to add the unpaired domains
            sum_seq_dist = accumulated_distance + abs(len(domIdxA) - len(domIdxB))
            # update the DDS or DDS_anchor and total domain counts depending on whether or not the domain is an anchor domain
            if shared_domain.split(".")[0] in anchor_domains:
                S_anchor += max(len(domIdxA),len(domIdxB))
                domain_difference_anchor += sum_seq_dist
            else:
                S += max(len(domIdxA),len(domIdxB))
                domain_difference += sum_seq_dist

        if S_anchor != 0 and S != 0:
            DDS_non_anchor = domain_difference / float(S)
            DDS_anchor = domain_difference_anchor / float(S_anchor)

            # Calculate proper, proportional weight to each kind of domain
            non_anchor_prct = S / float(S + S_anchor)
            anchor_prct = S_anchor / float(S + S_anchor)

            # boost anchor subcomponent and re-normalize
            non_anchor_weight = non_anchor_prct / (anchor_prct * anchorweight + non_anchor_prct)
            anchor_weight = anchor_prct * anchorweight / (anchor_prct * anchorweight + non_anchor_prct)

            # Use anchorweight parameter to boost percieved rDDS_anchor
            DDS = (non_anchor_weight * DDS_non_anchor) + (anchor_weight * DDS_anchor)

        elif S_anchor == 0:
            DDS_non_anchor = domain_difference / float(S)
            DDS_anchor = 0.0
            DDS = DDS_non_anchor

        else:  # only anchor domains were found
            DDS_non_anchor = 0.0
            DDS_anchor = domain_difference_anchor / float(S_anchor)
            DDS = DDS_anchor

        DDS = 1 - DDS  # transform into similarity
        # GK INDEX
        #  calculate the Goodman-Kruskal gamma index
        Ar = [item for item in A_list]
        Ar.reverse()
        GK = max([calculate_GK(A_list, B_list, nbhood), calculate_GK(Ar, B_list, nbhood)])

        Distance = 1 - (Jaccardw * Jaccard) - (DDSw * DDS) - (GKw * GK)
        if Distance < 0:
            print("Negative distance detected!")
            print("J: " + str(Jaccard) + "\tDDS: " + str(DDS) + "\tGK: " + str(GK))
            print("Jw: " + str(Jaccardw) + "\tDDSw: " + str(DDSw) + "\tGKw: " + str(GKw))
            sys.exit()
    except KeyError:
        print A, B,shared_domain
    return Distance, Jaccard, DDS, GK, DDS_non_anchor, DDS_anchor, S, S_anchor