#scale up faction 
#x=0.8 for membrane technologies 
#x=0.6 for volume increase technologies  

def scaleup_eq(C1,M1,M2,tec):
    if tec=="NF" or tec=="EDBM" or tec=="MED":
        x=0.6  
    else:
        x=0.6        
    C2=C1*(M2/M1)**x
    return C2

def scaleup(C1,M1,M2):
    x=0.6        
    C2=C1*(M2/M1)**x
    return C2