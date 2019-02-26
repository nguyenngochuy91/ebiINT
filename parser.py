#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Parse data from EBI INT, given a operon text file (Store operon gene name), given intact text file (protein 
              interaction file)
    Start   : 
    End     : 
'''

column = {4:"Alias(es) interactor A",
          5:"Alias(es) interactor B",6:'Interaction detection method(s)',
          9: 'Taxid interactor A',10:'Taxid interactor B',11: 'Interaction type(s)',
          14: 'Confidence value(s)'}
# this class take in intact.txt, and a species name and get all the protein interaction
def getGeneName(string):
    string = string.split("|")[1:]
    for item in string:
        item = item.split(":")
        dbName = item[0]
        name   = item[1]
        if dbName=="uniprotkb":
            if "(gene name)" in name:
                return name.replace("(gene name)","")
        
class MITAB(object):
    def __init__(self,intactFile,speciesName="ecoli"):
        self.dictionary  = {}
        self.process(speciesName,intactFile)
    def process(self,speciesName,intactFile): 
        handle = open(intactFile,"r")
        # get rid of first line
        handle.readline()
        line   = handle.readline().strip()
        count  = 0
        while line:
            item = line.split("\t")
            info = {}
            for i in range(len(item)):
                try:
                    info[column[i]]= item[i] # only store neccessary info
                except:
                    pass
            # take a look at taxid interactor A and B, if they are from the same speciesName
            if speciesName in info['Taxid interactor A'] and speciesName in info['Taxid interactor B']:
                # work on alias name to get the uniprotKB gene name
                nameA = getGeneName(info["Alias(es) interactor A"])
                nameB = getGeneName(info["Alias(es) interactor B"])
                info['Confidence value(s)'] = float(info['Confidence value(s)'].split(":")[-1])
                info['Interaction type(s)']  = info['Interaction type(s)'].split("(")[-1][:-1]
#             MITAB   if nameA=="atpA" and nameB=="atpD":
#                    print (info)
#                    break
                # generate a gene Name to store new info
                if nameA not in self.dictionary:
                    self.dictionary[nameA] = {}
                # if both interaction of gene A, gene B already store, pass, else, do this
                if nameB not in self.dictionary[nameA]:
                    self.dictionary[nameA][nameB]= {'Interaction detection method(s)':info['Interaction detection method(s)'],
                                           'Interaction type(s)':info['Interaction type(s)'],
                                           'Confidence value(s)':info['Confidence value(s)']}
                    count+=1
            line = handle.readline().strip()
            
        handle.close()
        print ("There are {} of protein interactions ".format(count))
    def getInteractionOfGene(self,geneName):
        try:
            return self.dictionary[geneName]
        except:
            pass
        
    def getInteractionOfGeneAGeneB(self,geneA,geneB):
        try:
            return self.dictionary[geneA][geneB]
        except:
            pass
    
            

                
#intact = MITAB("intact.txt")    