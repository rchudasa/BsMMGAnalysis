import json,sys,os


listOfDSets={
'ParkingBPH_3D':'/ParkingBPH3/Run2018D-20Jun2021_UL2018-v1/AOD',
}

for stat in listOfDSets:
    print(" Doing for "+stat+" : "+listOfDSets[stat] )

    os.system('set -x ; dasgoclient --query="run lumi dataset='+listOfDSets[stat]+'" > tmp ; set +x;')
    f=open('tmp','r')
    txt=f.readlines()
    f.close()
    strToPrint="{\n"
    for l in txt:
        strToPrint+="\n"
        l=l[:-1]
        run,lumilist=l.split(' ')
        run=run.replace(" ","")
        lumilist=lumilist.replace('[','')
        lumilist=lumilist.replace(']','')
        lumilist=lumilist.split(',')
        strToPrint+='"'+run+'" : ['
        for lumi in lumilist[:-1]:
            strToPrint+="["+lumi+","+lumi+"],"
        lumi=lumilist[-1]
        strToPrint+="["+lumi+","+lumi+"] ],"
    strToPrint=strToPrint[:-1]+"\n\n}"
    ofname=stat+'.json'
    f=open(ofname,'w')
    f.write(strToPrint)
    f.close()
    os.system('compareJSON.py --and '+ ofname + ' '+ ofname +' '+ ofname) 
