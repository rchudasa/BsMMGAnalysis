def getCustomization(bmmgTemplate,minDimuonMass,maxDimuonMass,maxDR,txtOffset=""):
    minDimuonMass_str="{0:0.2f}".format(minDimuonMass)
    maxDimuonMass_str="{0:0.2f}".format(maxDimuonMass)
    maxDR_str        ="{0:0.2f}".format(maxDR)
    tagname="BsToMMGIMin"+minDimuonMass_str.replace(".","p")
    tagname+=      "IMax"+maxDimuonMass_str.replace(".","p")
    tagname+=      "dRMax"+maxDR_str.replace(".","p")
#     print(tagname)
#     return
    modifiedCFF=""
    for l in bmmgTemplate:
        lm=l.replace("BsToMMG",tagname)
        if "MaxInvMass = cms.vdouble( 6.0 )" in l:
            lm=l.replace("6.0",maxDimuonMass_str)
        if "MinInvMass = cms.vdouble( 4.5 )" in l:
            lm=l.replace("4.5",minDimuonMass_str)
        if "MaxDr = cms.double( 2.0 )" in l:
            lm=l.replace("2.0",maxDR_str)
        modifiedCFF+=txtOffset+lm
    pathName="process.HLT_DoubleMu4_3_"+tagname+"_v0"
    return pathName,modifiedCFF

f=open("bs2mmgParts.py",'r')
bmmgParts=f.readlines()
f.close()

minDimuonMasses=[2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25]
maxDimuonMasses=[4.5,6.0]
maxDR=[0.75,1.0,1.25,1.5,1.75,2.0]
print("Generating ",len(maxDimuonMasses),"*",len(minDimuonMasses),"*",len(maxDR),\
      " : ",len(minDimuonMasses)*len(maxDimuonMasses)*len(maxDR)," Paths")

hlt_cff ="import FWCore.ParameterSet.Config as cms\n\n"
hlt_cff+="def customizeForHLTdev_Bs2MMG(process):\n"
pname_list=""
for maxDM in maxDimuonMasses:
    for maxdr in maxDR:
        for minDM in minDimuonMasses:
            pname,cff_tmp=getCustomization(bmmgParts,minDM,maxDM,maxdr,"    ")
            hlt_cff+=cff_tmp
            pname_list+="    " + "process.HLTSchedule += cms.Schedule(*[ "+pname + " ] )\n"
#             print(pname)
print("making file  hltCustomizationForTriggerStudiesBMMG.py  ")
f=open("hltCustomizationForTriggerStudiesBMMG.py",'w')
f.write(hlt_cff)
f.write('\n\n')
f.write(pname_list)
f.write('\n'+"    "+"return process\n\n")
f.close()

