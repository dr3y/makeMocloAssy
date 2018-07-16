import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import os
import math
#import beeswarm as bs
import sys
import time
import pydna
import itertools as it
import datetime
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.assembly import Assembly as pydAssembly
from Bio.Restriction import BsaI
from Bio.Restriction import BbsI
from Bio.Restriction import AarI
from copy import deepcopy as dc
import ipywidgets as widgets
from collections import defaultdict

enzymes = \
          {"BsaI":BsaI,
           "BbsI":BbsI,
           "AarI":AarI}
enlist = ["BsaI",
          "BbsI",
          "AarI",
          "gibson"]
prevplate = None
selenzyme = "gibson" #which enzyme to assemble everything with
chewnt = 40
frags = [] #fragments in the reaction
#the following lists the components in each well, in uL
gga = \
[["component","volume"],
 #["buffer10x",0.4],
 #["ATP10mM",0.4],
 #["BsaI", 0.2],
 #["ligase",0.2],
 ["NEBbuffer",0.4],
 ["NEBenzyme",0.2],
 ["water",1.4],
 ["dnasln",1],
 ]
gibassy = \
[["component","volume"],
["GGAMM",1],
["dnasln",1]]
ggaPD = pd.DataFrame(gga[1:],columns=gga[0]) #this just turns it into a data frame
gibassyPD = pd.DataFrame(gibassy[1:],columns=gibassy[0])

ggaFm = 8.0
ggavecGm = 8.0
gibFm = 6.0
gibvecFm = 6.0
partsFm = ggaFm  #default is gga
vectorFm = ggavecGm
source = "384PP_AQ_BP"
ptypedict = {
            "ASSGGA04":"384PP_PLUS_AQ_BP",
            "ASSGIB01":"384LDV_PLUS_AQ_BP",
            "ASSGIB02":"384PP_AQ_BP"}
waterwell = "P23" #in your source plate, include one well that is just full of water.
dnaPath = ".\\DNA\\"

#go down and look at makeEchoFile

def startText():
    print("Welcome to Moclo Assembly Helper V1")
    print("===================================")

def pickEnzyme():
    """asks the user about what kind of enzyme s/he wants to use"""
    print("Which enzyme would you like to use?")
    for el in range(len(enlist)):
                print("[{}]  {}".format(el,enlist[el]))
    print()
    userpick = int(input("type the number of your favorite! "))
    selenzyme = enlist[userpick].lower()
    print("===================================")
    return selenzyme
def findExpts(path):
    """gets a list of files/folders present in a path"""
    walkr = os.walk(path)
    dirlist = [a for a in walkr]
    expts = []
    #print(dirlist)
    #for folder in dirlist[1:]:
    folder = ['.']
    for fle in dirlist[0][2]:
        if(fle[-3:]=='csv'):
            try:
                fline = open(folder[0]+'\\'+fle,'r').readline().split(',')
                if("promoter" in fline):
                    expts+=[('{}\\{}'.format(folder[0],fle),fle[:-4])]
            except IOError:
                pass
        if(fle[-4:]=='xlsx'):
            try:
                xl_file = pd.ExcelFile(folder[0]+'\\'+fle)
                dfs = {sheet_name: xl_file.parse(sheet_name)
                          for sheet_name in xl_file.sheet_names}
                #print(dfs.keys()
                if(dfs["Sheet1"].columns[0] == "promoter"):
                    expts+=[('{}\\{}'.format(folder[0],fle),fle[:-5])]
            except (IOError,KeyError) as e:
                pass
    return sorted(expts)[::-1]

def findPartsLists(path):
    """gets a list of files/folders present in a path"""
    walkr = os.walk(path)
    dirlist = [a for a in walkr]
    #print dirlist
    expts = []
    for fle in dirlist[0][2]:
        #print fle
        if(fle[-4:]=='xlsx'):
            try:
                xl_file = pd.ExcelFile(path+'\\'+fle)
                dfs = {sheet_name: xl_file.parse(sheet_name)
                          for sheet_name in xl_file.sheet_names}
                #print(dfs.keys()
                if("parts" in list(dfs.keys())[0]):
                    expts+=[('{}\\{}'.format(path,fle),fle[:-4])]
            except IOError:
                pass
    return sorted(expts)[::-1]

def pickPartsList():
    """user interface for picking a list of parts to use. This list must
    contain the concentration of each part as well as the 384 well location
    of each part at minimum, but better to have more stuff. Check my example
    file."""
    print("Searching for compatible parts lists...")
    pllist = findPartsLists(".\\partslist")
    pickedlist = ''
    if(len(pllist) <=0):
        print("could not find any parts lists :(. Make sure they are in a \
                seperate folder called 'partslist' in the same directory as this script")
    else:
        print("OK! I found")
        print()
        for el in range(len(pllist)):
            print("[{}]  {}".format(el,pllist[el][1]))
        print()
        if(len(pllist)==1):
            pickedlist = pllist[0][0]
            print("picked the only one in the list!")
        else:
            userpick = int(input("type the number of your favorite! "))
            pickedlist = pllist[userpick][0]
    openlist = pd.ExcelFile(pickedlist)
    print("===================================")
    return openlist

def pickAssembly():
    """user interface for defining assemblies to build"""
    #manual = raw_input("would you like to manually enter the parts to assemble? (y/n)")

    manual = "n"
    if(manual == "n"):
        print("searching for compatible input files...")
        time.sleep(1)
        pllist = findExpts(".")
        #print pllist
        pickedlist = ''
        if(len(pllist) <=0):
            print("could not find any assembly files")
        else:
            print("OK! I found")
            print()
            for el in range(len(pllist)):
                print("[{}]  {}".format(el,pllist[el][1]))
            print()
            if(len(pllist)==1):
                pickedlist = pllist[0][0]
                print("picked the only one in the list!")
            else:
                userpick = int(input("type the number of your favorite! "))
                pickedlist = pllist[userpick][0]
        openlist = pd.read_csv(pickedlist)
        print("===================================")
        return openlist,pickedlist
    else:
        print("sorry I haven't implemented this yet")
        pickAssembly()
    return pd.read_csv(aslist),aslist
def echoline(swell,dwell,tvol,sptype = source,spname = "Source[1]",\
                                    dpname = "Destination[1]",platebc=""):
    #if(platebc!=""):
    #    sptype = ptypedict[platebc]
    return "{},{},{},{},,,,,{},{},{}\n".format(spname,platebc,sptype,swell,\
                                                            dpname,dwell,tvol)
def echoSinglePart(partDF,partname,partfm,dwell):
    """calculates how much of a single part to put in for a number of fm."""
    try:
        pwell = partDF[partDF.part==partname].well.iloc[0]
    except IndexError:
        raise ValueError("Couldn't find the right part named '"+\
          partname+"'! Are you sure you're using the right parts list?")
        return None, None, None
    pDseq = makeDseqFromDF(partname,partDF)
    pconc = partDF[partDF.part==partname]["conc (nM)"]
    #concentration of said part, in the source plate
    if(len(pconc)<=0):
        #in this case we could not find the part!
        raise ValueError("Part "+part+" had an invalid concentration!"+\
                            " Are you sure you're using the right parts list?")
    pconc = pconc.iloc[0]
    pplate = partDF[partDF.part==partname]["platebc"].iloc[0]
    platet = partDF[partDF.part==partname]["platetype"].iloc[0]
    e1,e2 = echoPipet(partfm,pconc,pwell,dwell,sourceplate=pplate,sptype=platet)
    return e1,e2,pDseq,pplate,platet
def echoPipet(partFm,partConc,sourcewell,destwell,sourceplate=None,sptype=None):
    """does the calculation to convert femtomoles to volumes, and returns the finished echo line"""
    pvol = (partFm/partConc)*1000
    evol = int(pvol)
    if(evol <= 25):#im not sure what happens when the echo would round to 0.
                    #better safe than sorry and put in one droplet.
        evol = 25
    if(sourceplate==None):
        print("===> transfer from {} to {}, {} nl".format(sourcewell,destwell,evol))
        echostring = echoline(sourcewell,destwell,evol)
    else:
        print("===> transfer from {}, plate {} to {}, {} nl".format(sourcewell,sourceplate,destwell,evol))
        echostring = echoline(sourcewell,destwell,evol,spname =sourceplate,sptype= sptype,platebc = sourceplate)
    return echostring, evol

def makeDseqFromDF(part,partslist,col = "part"):
    """looks up the part named "part" in the column specified as col, and converts it into a pydna object."""
    pseq = partslist[partslist[col] == part].sequence.iloc[0].lower()
    pcirc = partslist[partslist[col] == part].circular.iloc[0]
    p5pover = int(partslist[partslist[col] == part]["5pend"].iloc[0])
    p3pover = int(partslist[partslist[col] == part]["3pend"].iloc[0])

    povhg = int(p5pover)
    pseqRC = str(Dseq(pseq).rc()).lower()
    if(p5pover > 0):
        pseq = pseq[p5pover:]
    elif(p5pover<0):
        pseqRC = pseqRC[:p5pover]
    if(p3pover <0):
        pseq = pseq[:p3pover]
    elif(p3pover >0):
        pseqRC = pseqRC[p5pover:]
    pDseq = Dseq(pseq,pseqRC,ovhg=povhg)
    #this defines a dsdna linear sequence
    if(pcirc):
        #this makes the sequence circular, if we have to
        pDseq = pDseq.looped()
    return pDseq
def bluntLeft(DSseq):
    """returns true if the left hand side of DSseq is blunt"""
    if(DSseq.ovhg == 0):
        return True
    else:
        return False
def bluntRight(DSseq):
    """returns true if the left hand side of DSseq is blunt"""
    top = len(DSseq.watson)
    bot = len(DSseq.crick)
    diff = top-bot+DSseq.ovhg
    if(diff == 0):
        return True
    else:
        return False
def isNewDseq(newpart,partlist):
    new = True
    dsnewpart = Dseqrecord(newpart)
    for part in partlist:
        if(len(part) != len(newpart)):
            continue
        dspart = Dseqrecord(part)
        if(dspart.linear==False):
            try:
                dspart = dspart.synced(dsnewpart)
            except Exception:
                continue
        if(dspart == dsnewpart):
            new=False
            break
    return new
def allCombDseq(partslist,resultlist = []):
    '''recursively finds all possible paths through the partslist'''
    if(len(partslist)==1):
        return partslist#[[a] for a in partslist[0]]
    else:
        result = []
        for p in range(len(partslist)):
            newplist = dc(partslist)
            part = newplist.pop(p)
            prevresult = allCombDseq(newplist)
            partstoadd = []
            for z in prevresult:
                if(isNewDseq(z,result)):
                    partstoadd+=[z]
            for b in prevresult:
                #try to join the given part to everything else
                if((not bluntRight(part)) and (not bluntLeft(b))):
                    #this means we don't allow blunt ligations!
                    try:
                        newpart= part+b
                        try:
                            if((not bluntRight(newpart)) and (not bluntLeft(newpart))):
                                newpart = newpart.looped()
                                #this thing will return TypeError if it can't be
                                #looped
                            #result+= [loopart]
                        except TypeError:
                            #this happens if the part can't be circularized
                            pass
                        if(isNewDseq(newpart,result)):
                            result+=[newpart]
                    except TypeError:
                        #this happens if the parts don't have the right sticky ends
                        pass
            result+=partstoadd
        return result

def pushDict(Dic,key,value):
    """adds a value to a dictionary, whether it has a key or not"""
    try:
        pval = Dic[key]
    except KeyError:
        if(type(value)==list or type(value)==tuple):
            value = tuple(value)
            pval = ()
        elif(type(value)==str):
            pval = ""
        elif(type(value)==int):
            pval = 0
        elif(type(value)==float):
            pval = 0.0
    Dic[key] =pval + value
def findFilesDict(path=".",teststr = "promoter"):
    """gets a list of files/folders present in a path"""
    walkr = os.walk(path)
    dirlist = [a for a in walkr]
    expts = {}
    #print(dirlist)
    #for folder in dirlist[1:]:
    folder = [path]
    #print(dirlist)
    for fle in dirlist[0][2]:
        if(fle[-3:]=='csv'):
            try:
                print('{}\\{}'.format(folder[0],fle))
                fline = open(folder[0]+'\\'+fle,'r').readline().split(',')
                if(teststr in fline):
                    expts[fle[:-4]]='{}\\{}'.format(folder[0],fle)
            except IOError:
                pass
        if(fle[-4:]=='xlsx'):
            try:
                xl_file = pd.read_excel(folder[0]+'\\'+fle)
                #dfs = {sheet_name: xl_file.parse(sheet_name)
                #          for sheet_name in xl_file.sheet_names}
                #print(dfs.keys()
                #print(xl_file.columns)
                if(teststr in xl_file.columns):
                    #print("found")
                    expts[fle[:-5]]='{}\\{}'.format(folder[0],fle)
            except (IOError,KeyError) as e:
                pass
    return expts
def findPartsListsDict(path,teststr = "parts_1"):
    """gets a list of files/folders present in a path"""
    walkr = os.walk(path)
    dirlist = [a for a in walkr]
    #print(dirlist[0][2])
    expts = {}
    for fle in dirlist[0][2]:
        #print fle
        if((fle[-4:]=='xlsx') or (fle[-4:]=='xlsm')):
            try:
                xl_file = pd.ExcelFile(path+'\\'+fle)
                dfs = {sheet_name: xl_file.parse(sheet_name)
                          for sheet_name in xl_file.sheet_names}
                #print(dfs)
                #print(dfs.keys())
                if(teststr in list(dfs.keys())[0]):
                    expts[fle[:-5]] = '{}\\{}'.format(path,fle)
            except IOError:
                pass
    return expts

def findDNAPaths(startNode,nodeDict,edgeDict):
    """given a start, a dictionary of nodes, and a dictionary of edges,
    find all complete paths for a DNA molecule
    Complete is defined as: producing a molecule with all blunt edges,
    or producing a circular molecule."""
    #we assemble the DNA sequences from left to right.
    nnode = dc(nodeDict)
    noderight = nnode[startNode][1] #the right-hand overhang of the node in question.
    del nnode[startNode]
    destinations = edgeDict[noderight] #this could contain only one entry, the starting node
    seqs = [] #haven't found any complete paths yet
    nopaths = True
    candidateSeqs = []
    if(noderight != "blunt"): #blunt cannot go on
        for destination in destinations:
            #go through the list of destinations and see if we can go forward
            if(destination[1]==0): #this node links to something else
                if(destination[0] in nnode): #we havent visited it yet
                    nopaths = False

                    newpaths = findDNAPaths(destination[0],nnode,edgeDict) #find all paths from there!
                    for path in newpaths:
                        candidateSeqs+=[[startNode]+path]
    if(nopaths): #if we dont find any paths, call it good
        candidateSeqs+=[[startNode]]
    #print("canseqs is {}".format(candidateSeqs))
    return candidateSeqs
def getOverhang(Dnaseq,side="left"):
    """extracts the overhang in the DNA sequence, either on the left or right sides.
    If the dna sequence is blunt, then the returned overhang is called 'blunt'"""

def DPallCombDseq(partslist):
    '''Finds all paths through the partsist using a graph type of approach.
    First a graph is constructed from all possible overhang interactions,
    then the program makes paths from every part to a logical conclusion
    in the graph, then it backtracks and actually assembles the DNA.'''
    #actually, we need to produce a graph which describes the parts FIRST
    #then, starting from any part, traverse the graph in every possible path and store
    #the paths which are "valid" i.e., produce blunt ended or circular products.
    edgeDict = defaultdict(lambda : []) #dictionary of all edges in the partslist!
    nodeDict = {}#defaultdict(lambda : [])
    partDict = {}#defaultdict(lambda : [])
    pind = 0
    for part in partslist:
        Lend = ""
        Rend = ""
        Ltype,Lseq = part.five_prime_end()
        Rtype,Rseq = part.three_prime_end()
        if(Ltype == "blunt"):
            Lend = "blunt"
            edgeDict[Lend].append([pind,0])
            #pushDict(edgeDict,Lend,((pind,0),))
        else:
            if(Ltype == "3'"):
                Lend = str(Dseq(Lseq).rc()).lower()
            else:
                Lend = str(Lseq).lower()
            edgeDict[Lend].append([pind,0])
            #pushDict(edgeDict,Lend,((pind,0),))
        if(Rtype == "blunt"):
            Rend = "blunt"
            edgeDict[Rend].append([pind,1])
            #pushDict(edgeDict,Rend,((pind,1),))
        else:
            if(Rtype == "5'"):
                Rend = str(Dseq(Rseq).rc()).lower()
            else:
                Rend = str(Rseq).lower()
            edgeDict[Rend].append([pind,1])
            #pushDict(edgeDict,Rend,((pind,1),))
        nodeDict[pind] = (Lend,Rend)
        pind+=1
    paths = []
    for pind in list(nodeDict.keys()):
        paths += findDNAPaths(pind,nodeDict,edgeDict)
    goodpaths = []
    #print("paths are {}".format(paths))
    for path in paths:
        #print("path is")
        #print(path)
        fpart = path[0]
        rpart = path[-1]
        npart = False
        if(nodeDict[fpart][0]=="blunt" and nodeDict[rpart][1]=="blunt"):
            #this means we have a blunt ended path! good
            npart = True
            accpart = partslist[fpart]
            for pind in path[1:]:
                accpart+=partslist[pind]

        elif(nodeDict[fpart][0]==nodeDict[rpart][1]):
            npart = True
            #this means we have a circular part! also good!
            accpart = partslist[fpart]
            for pind in path[1:]:
                accpart+=partslist[pind]
            accpart=accpart.looped()
        if(npart):
            if(isNewDseq(accpart,goodpaths)):
                goodpaths+=[accpart]

    return goodpaths
def chewback(seqtochew,chewamt,end="fiveprime"):
    """chews back the amount mentioned, from the end mentioned."""
    wat = seqtochew.watson
    cri = seqtochew.crick

    if(len(seqtochew) > chewamt*2+1):
        if(end=="fiveprime"):
            cwat = wat[chewamt:]
            ccri = cri[chewamt:]

        else:
            cwat = wat[:-chewamt]
            ccri = cri[:-chewamt]
        newseq = Dseq(cwat,ccri,ovhg = chewamt)
        return newseq
    else:
        return None

def makeEchoFile(parts,aslist,gga=ggaPD,partsFm=partsFm,source=source,\
            output = "output.csv",selenzyme=selenzyme,fname="recentassembly",\
            protocolsDF=None,sepfiles=True,sepfilename="outputLDV.csv"):
    """makes an echo csv using the given list of assemblies and source plate of
    parts..
    inputs:
        parts: dataframe of what's in the source plate
        aslist: dataframe of what we need to assemble
        gga: a short dictionary indicating what volume of all the components
            go into the reaction mix
        partsFm: how many femtomoles of each part to use
        source: the name of the source plate. like "384PP_AQ_BP or something
        output: the name of the output file
        selenzyme: the enzyme we are going to use for assembly. everything
            is assembled with the same enzyme! actually this does nothing because
            the enzyme is taken from the aslist thing anyway
        fname: this is the name of the folder to save the successfully assembled
            dna files into
        protocolsDF: a dataframe containing a descriptor for different possible
            protocols. For instance it would say how much DNA volume and
            concentration we need for GGA or gibson."""

    #this is the boilerplate columns list
    outfile = "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,\
    Sample ID,Sample Name,Sample Group,Sample Comment,Destination Plate Name,\
    Destination Well,Transfer Volume\n"
    outfile2 = "Source Plate Name,Source Plate Barcode,Source Plate Type,Source Well,\
    Sample ID,Sample Name,Sample Group,Sample Comment,Destination Plate Name,\
    Destination Well,Transfer Volume\n"
    f2init = len(outfile2)
    #this iterates through rows in the assembly list file. Each row
    #defines an assembly, with the columns representing what parts go in.
    #this may not be ideal but it's fairly human readable and we only do
    #four parts + vector for each assembly.
    fname = fname.split("\\")[-1]
    if("." in fname):
        fname = fname[:fname.index(".")]
        print("saving in folder ./DNA/{}".format(fname))
    #the following is for making a spreadsheet style sequence list for performing further assemblies
    prodSeqSpread = "well,part,description,type,left,right,conc (nM),date,numvalue,sequence,circular,5pend,3pend,length\n"
    prevplate = None
    prevtype = None
    for assnum in range(len(aslist)):
        assembly = aslist[assnum:assnum+1] #cuts out one row of dataframe
        dwell = assembly.targwell[assembly.targwell.index[0]] #well where assembly will happen

        #print("pick enzyme")
        #print(assembly)
        try:
            selenzyme = assembly.enzyme[assembly.enzyme.index[0]]
            #if the user forgot to define an enzyme assume it is BsaI. That's the most common one we use
        except KeyError:
            selenzyme = "BsaI"
        if(protocolsDF!=None):
            cprt_temp = "gga"
            if(selenzyme == "gibson"):
                cprt_temp = "gibson"
            curprot = {"dnasln": protocolsDF[(protocolsDF.protocol==cprt_temp)&\
                            (protocolsDF.component == "dnasln")].amount.iloc[0]}
            partsFm = curprot[curprot.component==partfm].amount.iloc[0]
            vectorFm = curprot[curprot.component==vectorfm].amount.iloc[0]
        else:
            curprot = ggaPD
            partsFm = ggaFm
            vectorFm = ggavecGm
            if(selenzyme == "gibson"):
                curprot = gibassyPD
                partsFm = gibFm
                vectorFm = gibvecFm
        water = float(curprot[curprot.component=="dnasln"].volume)*1000 #total amount of water, to start with
        print("assembling with "+selenzyme)
        aind = assembly.index[0] #necessary for dataframes probably because I'm dumb
        frags = []
        if(not selenzyme == "gibson"):
            enzyme = enzymes[selenzyme]
            esite = enzyme.site.lower()
            esiterc = str(Dseq(enzyme.site).rc()).lower()
        for col in assembly:
            if(col=="targwell"):#since every row is terminated by the "target well",
                                #we'll take this opportunity to put in the water
                if(int(water) <25):
                    #echo gets mad if you tell it to pipet significantly less than 25 nl
                    water = 25
                ewat = int(water) #the echo automatically rounds to the nearest 25,
                                #so it's not really necessary to round here.
                #dsrfrags = [Dseqrecord(a) for a in frags]
                #x = pydAssembly(dsrfrags,limit = 4)
                #print(frags)
                #print(len(frags))
                allprod= []
                nefrags = []
                cutfrags = []
                if(selenzyme != "gibson"):
                    enzyme = enzymes[selenzyme]
                for frag in frags:
                    if(selenzyme == "gibson"):
                        if(len(frag)>chewnt*2+1):
                            nefrags += [chewback(frag,chewnt)]
                        else:
                            raise ValueError("part with sequence "+frag+" is too "+\
                                            "short for gibson! (<= 80 nt)")
                    else:
                        newpcs = frag.cut(enzyme)
                        if(len(newpcs) == 0):
                            newpcs+=[frag]
                        for pcs in newpcs:
                            if(pcs.find(esite)+pcs.find(esiterc)==-2):
                                nefrags+=[pcs]
                allprod = DPallCombDseq(nefrags)
                goodprod = []
                newpath = dnaPath+fname
                Cname = ""
                try:
                    #this part gathers the "name" column to create the output sequence
                    Cname = assembly.name[assembly.name.index[0]]
                except KeyError:
                    Cname = ""
                if(Cname == "" or str(Cname) == "nan"):
                    Cname = "well"+dwell
                print("Parts in construct {}".format(Cname))
                if not os.path.exists(newpath):
                    print("made dirs!")
                    os.makedirs(newpath)

                num = 0
                for prod in allprod:
                    Cnamenum = Cname
                    if(len(allprod) > 1):
                        wout = open(newpath+"\\"+Cname+"_"+str(num)+".gbk","w")
                        Cnamenum = Cname+"_"+str(num)
                    else:
                        wout = open(newpath+"\\"+Cname+".gbk","w")
                    if(bluntLeft(prod) and bluntRight(prod)):
                        num+=1
                        goodprod+=[prod]
                        topo = ["linear","circular"][int(prod.circular)]
                        booltopo = ["FALSE","TRUE"][int(prod.circular)]
                        #wout.write("\r\n>Construct"+str(num)+"_"+topo)
                        un_prod = "_".join(Cnamenum.split())
                        wout.write("LOCUS       {}                {} bp ds-DNA     {} SYN 01-JAN-0001\n".format(un_prod,len(prod),topo))
                        wout.write("ORIGIN\n")
                        wout.write(str(prod)+"\n//")
                        now = datetime.datetime.now()
                        nowdate = "{}/{}/{}".format(now.month,now.day,now.year)
                        prodSeqSpread += "{},{},assembled with {},,,,30,{},,{},{},{},{},{}\n".format(\
                                        dwell,un_prod,          selenzyme,nowdate,prod,booltopo,0,0,len(prod))
                    wout.close()
                assembend = ["y","ies"][int(len(goodprod)>1)]
                print("Detected {} possible assembl{}".format(len(goodprod),assembend))
                frags = []
                if(water <=0):
                    print("WARNING!!!! water <=0 in well {}".format(dwell))
                else:
                    #print("water from {} to {}, {} nl".format(waterwell,dwell,ewat))
                    if(prevplate == None):
                        #print("normalwater")
                        outfile += echoline(waterwell,dwell,ewat)
                    else:
                        #print("platewater")
                        outfile += echoline(waterwell,dwell,ewat,spname =prevplate,sptype=prevtype,platebc = prevplate)
                    #add water to the well!
                print("")
            elif(col == "comment"):#skip this column!
                pass
            elif(col == "enzyme"):
                pass
            elif(col == "name"):
                pass
            else:
                part = assembly[col][aind]

                #print(assembly[col][aind])
                #print("currently on "+part)
                #part = assembly[col][aind] #well corresponding to the part we want
                if(str(part) == 'nan'):
                    #this means we skip this part
                    print("skip one!")
                else:
                    part = assembly[col][aind]
                    #this is the name of the part!
                    #parts[parts.part==assembly[col][aind]].well.iloc[0]
                    evol = 0
                    if(':' in str(part)):
                        #this means we have multiple parts to mix!
                        subparts = part.split(':')
                        t_partsFm = partsFm/len(subparts)
                        t_vecFm = vectorFm/len(subparts)
                        for subpart in subparts:
                            useFm = t_partsFm
                            if(col == "vector"):
                                #use the vector at lower concentration!!
                                useFm = t_vecFm
                            e1,e2,pDseq,prevplate,prevtype = echoSinglePart(parts,\
                                                        subpart,useFm,dwell)
                            frags+=[pDseq]
                            evol += e2
                            if(sepfiles):
                                if("LDV" in e1):
                                    outfile2+=e1
                                else:
                                    outfile+= e1
                            else:
                                outfile+= e1


                    else:
                        useFm = partsFm
                        if(col == "vector"):
                            #use the vector at lower concentration!!
                            useFm = vectorFm
                        e1,e2,pDseq,prevplate,prevtype = echoSinglePart(parts,\
                                                    part,useFm,dwell)
                        frags+=[pDseq]
                        evol += e2
                        if(sepfiles):
                            if("LDV" in e1):
                                outfile2+=e1
                            else:
                                outfile+= e1
                        else:
                            outfile+= e1
                    water=water-evol
    pspread = open(newpath+"\\"+fname+".csv","w")
    pspread.write(prodSeqSpread)
    pspread.close()
    seqdispDF = pd.read_csv(newpath+"\\"+fname+".csv",usecols=["well","part","circular","length"])
    display(seqdispDF)
    ofle = open(output,"w")
    ofle.write(outfile)
    ofle.close()
    if(sepfiles and (len(outfile2) > f2init)):
        print("wrote LDV steps in {}".format(sepfilename))
        ofle2 = open(sepfilename,"w")
        ofle2.write(outfile2)
        ofle2.close()
outitems = []


class assemblyFileMaker():
    def __init__(self,mypath="."):
        self.parts = findPartsListsDict(mypath+"\\partslist")
        #txtdisabl = False
        assemblies = []
        self.fname1 = widgets.Text(
            value="untitled",
            placeholder = "type something",
            description='Assembly File Name:',
            disabled=False
        )
        self.drop2 = widgets.Dropdown(
            options=self.parts,
            width=100,
            #value=2,
            description='parts list:',
        )
        self.but = widgets.Button(
            description='Select',
            disabled=False,
            button_style='', # 'success', 'info', 'warning', 'danger' or ''
            tooltip='Click me',

            #icon='check'
        )
        self.but.on_click(self.on_button_clicked)
        self.cbox = widgets.HBox([self.fname1,self.drop2,self.but])
        display(self.cbox)
    def add_row(self,b):
        self.addWidgetRow(labonly=False)
        outcols = [widgets.VBox(a) for a in self.outitems ]
        self.bigSheet.children=outcols
        b.disabled=True
        #print(b)
    def remove_row(self,b):
        print(b)
    def on_button_clicked(self,b):
        #txtdisabl = True
        b.disabled=True
        xl_file = pd.ExcelFile(self.drop2.value)
        dfs = {sheet_name: xl_file.parse(sheet_name)
                          for sheet_name in xl_file.sheet_names}
        sheetlist = list(dfs.keys())
        self.p = pd.DataFrame.append(dfs["parts_1"],dfs["Gibson"])
        self.collabels = ["Promoter","UTR","CDS","Terminator","vector1","vector2","name",""]
        self.outitems = []
        self.ddlay = widgets.Layout(width='75px',height='30px')
        self.eblay = widgets.Layout(width='50px',height='30px')
        self.currow = 0
        self.addWidgetRow(labonly=True)
        self.addWidgetRow(labonly=False)
        #outitems+=[[labwidg,interwidg]]
        outcols = [widgets.VBox(a) for a in self.outitems ]
        self.bigSheet=widgets.HBox(outcols)
        display(self.bigSheet)
        self.fname1.disabled = True
        #self.outcolnum = 0
    def addWidgetRow(self,labonly=True):
        self.currow +=1
        outcolnum=0
        for col in self.collabels:
            if(labonly):
                interwidg = widgets.Label(col)
            else:

                if(col=="name"):
                    interwidg = widgets.Text(\
                            layout=self.ddlay)
                elif(col==""):
                    but1 = widgets.Button(\
                        description='+',
                        button_style='success',
                        tooltip='row '+str(self.currow),
                        layout=self.eblay
                    )
                    but2 = widgets.Button(\
                        description='-',
                        button_style='danger',
                        tooltip='row '+str(self.currow),
                        layout=self.eblay
                    )
                    but1.on_click(self.add_row)
                    but2.on_click(self.remove_row)
                    interwidg =widgets.HBox([but1,but2])
                else:
                    if("vector" in col):
                        coltext = "UNS"
                    else:
                        coltext = col
                    interwidg = widgets.Dropdown(\
                            options=list(self.p[self.p.type==coltext].part)+["None"],\
                            layout=self.ddlay)
            try:
                self.outitems[outcolnum]+=[interwidg]
            except IndexError:
                self.outitems+=[[interwidg]]
            outcolnum +=1


def makeAssemblyFile(mypath="."):
    """this function will assist the user with making assembly .csv files!"""
    x=assemblyFileMaker(mypath=".")

def makeInteractive(mypath="."):
    oplist = findFilesDict(mypath+"\\assemblies")
    parts = findPartsListsDict(mypath+"\\partslist")

    drop1 = widgets.Dropdown(
        options=oplist,
        #value=2,
        description='Assembly:',
    )
    drop2 = widgets.Dropdown(
        options=parts,
        #value=2,
        description='parts list:',
    )
    but = widgets.Button(
        description='Select',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        #icon='check'
    )

    #button = widgets.Button(description="Click Me!")
    #display(button)
    #print(oplist)
    def on_button_clicked(b):
        xl_file = pd.ExcelFile(drop2.value)
        #print(drop1.value)
        if(drop1.value[-4:]=="xlsx" or drop1.value[-3:]=="xls"):
            x=pd.read_excel(drop1.value)
        else:
            x=pd.read_csv(drop1.value)
        dfs = {sheet_name: xl_file.parse(sheet_name)
                          for sheet_name in xl_file.sheet_names}
        sheetlist = list(dfs.keys())
        p = pd.DataFrame.append(dfs["parts_1"],dfs["Gibson"])

        makeEchoFile(p,x,fname = drop1.value, \
                    output = ".\\output\\output.csv",\
                    sepfilename=".\\output\\outputLDV.csv")

        #print(drop1.value+" and "+drop2.value)

    but.on_click(on_button_clicked)
    cbox = widgets.HBox([drop1,drop2,but])
    display(cbox)


def runProgram():
    #x=pd.read_csv(insheet,sep=",")
    #pickhand = raw_input("is this for the echo? (y/n)")
    pickhand = 'y'
    xl_file=pickPartsList()
    x,fname=pickAssembly()
    #enz=pickEnzyme()
    #p=pd.read_csv("partslist/CIDAR_parts_plate_ASS.csv",sep=",")


    #pd.ExcelFile("partslist/CIDAR_parts_plate_ASS.xlsx")
    dfs = {sheet_name: xl_file.parse(sheet_name)
          for sheet_name in xl_file.sheet_names}
    sheetlist = list(dfs.keys())
    p = pd.DataFrame.append(dfs["parts_1"],dfs["Gibson"])
    #print(p)
    try:
        if(pickhand=="n"):
            makeHandFile(p,x)
        else:
            makeEchoFile(p,x,fname = drop1.value, \
                        output = ".\\output\\output.csv",\
                        sepfilename=".\\output\\outputLDV.csv")
    except ValueError as error:
        print("=========ERROR========")
        print(error)
        print ("")
        return False
    return True

if(__name__=="__main__"):
    progran = False
    while(not progran):
        progran = runProgram()
        if( not progran):
            z=input("Press Enter to try again")
        print("======================")
        print("")
        print("")
    z=input("Done! Press Enter to exit")




#def askForParts():
#    print "Searching for compatible parts lists..."
