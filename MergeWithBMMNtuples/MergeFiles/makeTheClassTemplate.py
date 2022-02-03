#!/usr/bin/env python3

fillingDebugOn =False
shortExt={\
'string'          :   'C',\
'Char_t'          :   'B',\
'UChar_t'         :   'b',\
'Short_t'         :   'S',\
'UShort_t'        :   's',\
'Int_t'           :   'I',\
'UInt_t'          :   'i',\
'Float_t'         :   'F',\
'Float16_t'       :   'f',\
'Double_t'        :   'D',\
'Double32_t'      :   'd',\
'Long64_t'        :   'L',\
'ULong64_t'       :   'l',\
'Long_t'          :   'G',\
'ULong_t'         :   'g',\
'Bool_t'          :   'O'\
}

maxLen='MAX_ARRAY_SIZE'

f=open('branchList','r')
l=f.readline()
treeName="_compiledTree"

branchDefenitionTxt=""
dataMemberDefenitionCount=""
dataMemberDefenitionArray=""

stlVectorToArrayConversionFuntion_Count=""
stlVectorToArrayConversionFuntion_Fill=""

fillArray_dict={}

while l:
    isAdress=False
    isArray =False
    arrayExt=""
    if l[:-1]=="":
        l=f.readline()
        continue
    if l[0]=="#":
        l=f.readline()
        continue
    items=l[:-1].replace(' ','').split(",")
    dataType = items[1]
    className=items[0]
    branchName=items[2].split('[')[0].replace("*","")
    outputBranchName = className + "_" + branchName
    
    if "*" in items[2] or "[" in items[2] : isAdress= True
    if "[" in items[2] : isArray= True
    
    if len(items)<4 and  isArray:
        print("Sanity Check Fali !! : ",l[:-1])
    if len(items)>3 : 
        arrayLengthName=className+"_"+items[3].replace('[','').replace(']','').strip()
       
    
    isVector=False
    if 'vector' in l:
        isVector=True
        baseDataType =dataType.replace("std::","").replace("vector","").replace("<","").replace(">","")
        baseDataType= baseDataType.replace("\t","")

        if "float"  in baseDataType :baseDataType='Float_t'
        if "double" in baseDataType :baseDataType='Double_t'
        if "int"    in baseDataType :baseDataType='Int_t'
        if "bool"   in baseDataType :baseDataType='Bool_t'
        
        
        if len(items)>3:
            if arrayLengthName not in fillArray_dict : fillArray_dict[arrayLengthName]={'slaveArrays':[],'sourceVects':[],'sizeFrom':""}
            fillArray_dict[arrayLengthName]['slaveArrays'].append(outputBranchName)
            if  fillArray_dict[arrayLengthName]['sizeFrom']=="":
                fillArray_dict[arrayLengthName]['sizeFrom']=className+"."+items[3].replace('[','').replace(']','')

        else:
            arrayLengthName=className+"_n"+branchName;
            if arrayLengthName not in fillArray_dict : fillArray_dict[arrayLengthName]={'slaveArrays':[],'sourceVects':[],'sizeFrom':""}
            
            tmp ="\tUInt_t "+arrayLengthName+";\n"
            dataMemberDefenitionArray+=tmp;

            tmp=""
            tmp+="\t"+arrayLengthName+" = "+className+"."+branchName+"->size() ;\n"
            tmp+="\t"+arrayLengthName+" = "+arrayLengthName+" > "+str(maxLen)+" ? "+str(maxLen)+" : "+arrayLengthName+" ;\n"
            stlVectorToArrayConversionFuntion_Count+=tmp
            fillArray_dict[arrayLengthName]['slaveArrays'].append(outputBranchName)
            fillArray_dict[arrayLengthName]['sizeFrom']=arrayLengthName

            branchDefenitionTxt+="\t"+treeName+'->Branch("'+arrayLengthName+'", &' +arrayLengthName+");\n"
        
        fillArray_dict[arrayLengthName]['sourceVects'].append(className+'.'+branchName)


        branchArray=outputBranchName

        tmp="\t" +baseDataType+"  "
        tmp+=branchArray+"["+str(maxLen)+"] ;\n"
        dataMemberDefenitionArray+=tmp;
        
       
        isAdress=True
        isArray =True
        dataType=baseDataType
        
    

    text="\t"+treeName+'->Branch("'+outputBranchName+'",'
    if  not isAdress : text+="&"
    if isVector:
        text+="\t( " + outputBranchName +" )"
    else:    
        text+="\t( " + className +"."+ branchName +" )"
    if isArray:
        arrayExt  =',"'+outputBranchName+"["
        arrayExt +=arrayLengthName
        arrayExt +="]"
        arrayExt +="/"+shortExt[dataType]+'"'

    text+=arrayExt
    text+=');'
    branchDefenitionTxt+=text+"\n"
#     print(text)
        
    l=f.readline()

# Fill the arrys here
stlVectorToArrayConversionFuntion_Fill+="\n\n"

for arryLens in fillArray_dict:
    tmp=""
    if(fillingDebugOn): tmp+='\tstd::cout<<__LINE__<<"'+arryLens+'"<<endl ; \n'
    tmp +="\tfor( UInt_t i=0; i< UInt_t ( "+fillArray_dict[arryLens]['sizeFrom']+") ; i++ ) { \n" 
    tmp +="\t\t if( i>= MAX_ARRAY_SIZE ) {    break;}\n "
    for arrs,src in zip(fillArray_dict[arryLens]['slaveArrays'],fillArray_dict[arryLens]['sourceVects']):

#        if(fillingDebugOn): tmp+='\tstd::cout<<__LINE__<<" '+arryLens+' "<<"'+src+'"<<endl ; \n'

        tmp+="\t\t"+arrs+"[i] = "+src+"->at(i);\n"
    tmp+="\n\t}\n\n"
    stlVectorToArrayConversionFuntion_Fill+=tmp
     
f=open("DataMembers.h",'w')
f.write("//      DATAMEMBERS   \n")
f.write(dataMemberDefenitionCount)
f.write("\n\n")
f.write(dataMemberDefenitionArray)
f.write("//      END DATA MEMBERS   \n\n")
f.close()


f=open("MemberFunction.h",'w')
f.write("void genericTreeBranchSelector::FillTheArraysFromVectors()\n")
f.write("{\n")
f.write(stlVectorToArrayConversionFuntion_Count)
f.write(stlVectorToArrayConversionFuntion_Fill)
f.write("\n}")
f.write("\n\n")
f.write("void genericTreeBranchSelector::setCompiledTreeBranches() ")
f.write("\n{\n\n")
f.write(branchDefenitionTxt)
f.write("\n}\n")
f.close()    
