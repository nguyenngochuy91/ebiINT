#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Parse data from EBI INT, given a operon text file (Store operon gene name),generate network for our graph
    Start   : 
    End     : 
'''
from parser import MITAB
import pydot
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
def parseOperonFile(operonFile):
    dictionary = {}
    handle     = open(operonFile,"r")
    for line in handle.readlines():
        line = line.strip().split()
        name = line[0]
        genes = line[1:]
        dictionary[name] = genes
    return dictionary

# generate graph, vertice are gene, edge are interaction, weight is the confidence value                    
def getOperonInteraction(operonDictionary,MITTABobject):
    scores = np.array([])
    operonName =[]
    for operon in operonDictionary:
        geneSet = operonDictionary[operon]
        operonName.append(operon)
        print ("Operon: {} with gene set {}".format(operon,geneSet))
        n = len(geneSet)
        total = n*(n-1)/2.0
        edgeWeights = 0.0
        # generate network
        graph = pydot.Dot(graph_type="graph")
        for i in range(len(geneSet)):
            newNode = pydot.Node(geneSet[i])
            graph.add_node(newNode)
            for j in range(len(geneSet)):
                if i!=j:
                    geneA = geneSet[i]
                    geneB = geneSet[j]
                    interaction = MITTABobject.getInteractionOfGeneAGeneB(geneA,geneB)
                    if interaction:
                        print ("Gene {} interact with gene {} as follow \n {}".format(
                                geneA,geneB,interaction))
                        oldE      = MITTABobject.getInteractionOfGeneAGeneB(geneB,geneA)
                        if oldE:
                            if oldE['Confidence value(s)']<interaction['Confidence value(s)']:
                                oldE[0].set('label',label ="weight={}\ntype={}".format(interaction['Confidence value(s)'],interaction['Interaction type(s)']))
                                edgeWeights= edgeWeights-float(oldE['Confidence value(s)'])+float(interaction['Confidence value(s)'])
                        else:
                            newE      = pydot.Edge(geneA,geneB,color="blue",label ="w={}\nt={}".format(interaction['Confidence value(s)'],interaction['Interaction type(s)']))
                            graph.add_edge(newE)
                            edgeWeights += float(interaction['Confidence value(s)'])
        score = edgeWeights/total
        scores = np.append(scores,score)
        graph.write_png("operon_graph/"+operon)
    # generate the z score
    zscores = st.zscore(scores)
    dictionary = {}
    for i in range(len(operonName)):
        dictionary[operonName[i]] = scores[i]
    return dictionary
# given the graph, compute a normalize score    

def parseNormalizeOperonFile(operonFile):
    name = []
    scores =[]
    dictionary = {}
    handle = open(operonFile,"r")
    handle.readline()
    for line in handle.readlines():
        item = line.strip().split()
        name.append(item[0])
        total      = float(item[1])
        scores.append(total)
    zscores = st.zscore(scores)
    dictionary = {}
    for i in range(len(name)):
        dictionary[name[i]] = scores[i]
    return dictionary
# given graphDictionary,eventDictionary, do analysis, draw plot and blah blah blah :)
def analysis(graphDictionary,eventDictionary,alpha):
    distance_graph = []
    for key in graphDictionary:
        distance_graph.append((eventDictionary[key],graphDictionary[key]))
    distance_graph.sort()
    x = [item[0] for item in distance_graph]
    y = [item[1] for item in distance_graph]

    # do the testing
    # pearson testing, assumption is that variables follow normal distributions
    # draw the 2 bar
    plt.bar(x,y,.1)
    plt.scatter(x,y)
    plt.plot(x,y)
    a = np.concatenate((x, y))
    k2, p = st.normaltest(a)
    print("p = {:g}".format(p)) 
    print ("null hypothesis: our data comes from a normal distribution, with threshold of {}".format(alpha))
    if p<alpha:
        print("The null hypothesis can be rejected")
    else:
        print("The null hypothesis cannot be rejected")
    corr,pval      = st.pearsonr(x,y)
    print ("Correlation of Pearson: {} with p : {}".format(corr,pval))
    rho, pvalSpear = st.spearmanr(x,y)
    print ("Correlation of Spearman: {} with p : {}".format(rho,pvalSpear))
    tau, pvalKen = st.kendalltau(x, y)
    print ("Correlation of Kendall: {} with p : {}".format(tau,pvalKen))
    
def main():
    print ("Parsing data form IntAct ...")
    intact = MITAB("intact.txt")
    print ("Generating graph, and provide normalize score of the graph interaction ...")
    graphDictionary= getOperonInteraction(parseOperonFile("gene_block_names_and_genes.txt"),intact)
    print (graphDictionary)
    print ("Parsing normalizing event score of each operon ...")
    eventDictionary = parseNormalizeOperonFile("operon_conservation_score/eColiConservedOperonsSorted.txt")
    # do pearson testing \
    analysis(graphDictionary,eventDictionary,1e-23)
#graphDictionary,eventDictionary = main()